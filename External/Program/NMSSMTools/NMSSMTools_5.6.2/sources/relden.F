      SUBROUTINE RELDEN(PAR,PROB)

**********************************************************************
*   Subroutine for the computation of the dark matter relic density
*   PROB(26) =/= 0  lightest neutralino is not LSP
*   PROB(30) =/= 0 excluded by WMAP (checked only if OMGFLAG=/=0)
*   PROB(30)  = -1 LSP is not the lightest neutralino in micrOMEGAs
*   PROB(30)  = -2 Problem in micrOMEGAs
*   PROB(31) =/= 0  excluded by DM SI WIMP-nucleon xs (checked if |OMGFLAG|=2 or 4)
*   PROB(61) =/= 0  excluded by DM SD WIMP-neutron xs (checked if |OMGFLAG|=2 or 4)
*   PROB(62) =/= 0  excluded by DM SD WIMP-proton xs (checked if |OMGFLAG|=2 or 4)
*   GRFLAG=0: the thermal relic density of the lightest neutralino is computed
*             using micrOMEGAs. Default value for NMHDECAY and NMSPEC
*   GRFLAG=1: if the gravitino is the LSP, first the thermal relic density of the
*             NLSP is computed using micrOMEGAs, then the gravitino non thermal
*             relic density is estimated as OMG(LSP) = OMG(NLSP)*mass(LSP)/mass(NLSP).
*             Default value in NMGMSB.
*
**********************************************************************

      IMPLICIT NONE

#include "micromegas/include/micromegas.fh"

      CHARACTER name*10,mess*20
      CHARACTER*8 process(5)

      INTEGER NORD(5),HORD(3),NBIN,OMGFLAG,MAFLAG,MOFLAG,GRFLAG
      INTEGER err,i,j,WW,fast,ok,pdg(5)

      DOUBLE PRECISION PAR(*),PROB(*)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW,PI
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION SMASS(3),AMASS(2),CMASS,SCOMP(3,3),PCOMP(2,2)
      DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION SST,SSB,SSL,COSB,SINB,TANB
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION tab(250),OMG,OMGMIN,OMGMAX,Xf
      DOUBLE PRECISION sigmaV,x(100),dNdx(100),EMIN,LAM
      DOUBLE PRECISION XENON_SI,XENON_SDn,XENON_SDp,PandaX_SI
      DOUBLE PRECISION LUX_SI,LUX_SDn,LUX_SDp,PICO60_SDp
      DOUBLE PRECISION CRESST_SI,DarkSide50_SI,lim
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION Nmass,SCcoeff,DELMB,Beps
      DOUBLE PRECISION pA0(2),pA5(2),nA0(2),nA5(2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5),BRSS(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION weight,vcsll,vcsbb
      DOUBLE PRECISION M32,CGR,MPL
      DOUBLE PRECISION scal

      COMMON/DELMB/DELMB
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/HIGGSPEC/SMASS,SCOMP,AMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/LAM/LAM
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRSS,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/M32/M32,CGR,MPL,GRFLAG

      DATA NORD/1,2,4,3,5/
      DATA HORD/2,1,3/
    
*   Test on LSP

      IF(MIN(DABS(MGL),DABS(MCHA(1)),MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .       MST1,MSB1,MSL1,MSNT,DABS(MNEU(1))).LT.511d-6)THEN
       PROB(26)=-DDIM(1d0,MIN(DABS(MGL),DABS(MCHA(1)),MUR,MUL,MDR,
     .       MDL,MLR,MLL,MNL,MST1,MSB1,MSL1,MSNT,DABS(MNEU(1)))/511d-6)
      ELSEIF(M32.GT.DABS(MNEU(1)))THEN
       PROB(26)=DDIM(1d0,MIN(DABS(MGL),DABS(MCHA(1)),MUR,MUL,MDR,
     .          MDL,MLR,MLL,MNL,MST1,MSB1,MSL1,MSNT)/DABS(MNEU(1)))
      ENDIF
      IF(PROB(26).NE.0d0)OMG=-1d0
      IF (OMGFLAG.EQ.0 .OR. PROB(26).NE.0d0) RETURN

*   Input parameters:

      PI=4d0*DATAN(1d0)
      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      CALL assignValW('alfSMZ',ALSMZ)
      CALL assignValW('MbMb',MB)
      CALL assignValW('Mtp',MT)

      CALL assignValW('At',PAR(12))
      CALL assignValW('Ab',PAR(13))
      CALL assignValW('Al',PAR(14))

      CALL assignValW('Lambda',PAR(1))
      CALL assignValW('Kappa',PAR(2))
      CALL assignValW('aLambda',PAR(5))
      CALL assignValW('aKappa',PAR(6))
      CALL assignValW('tB',TANBQ)
      CALL assignValW('mu',PAR(4))

      CALL assignValW('MG1' ,PAR(20))         
      CALL assignValW('MG2' ,PAR(21))         
      CALL assignValW('MG3' ,PAR(22))         
      CALL assignValW('Ml2' ,DSQRT(PAR(18)))  
      CALL assignValW('Ml3' ,DSQRT(PAR(10)))  
      CALL assignValW('Mr2' ,DSQRT(PAR(19)))  
      CALL assignValW('Mr3' ,DSQRT(PAR(11)))  
      CALL assignValW('Mq2' ,DSQRT(PAR(15)))  
      CALL assignValW('Mq3' ,DSQRT(PAR(7)))   
      CALL assignValW('Mu2' ,DSQRT(PAR(16)))  
      CALL assignValW('Mu3' ,DSQRT(PAR(8)))   
      CALL assignValW('Md2' ,DSQRT(PAR(17)))  
      CALL assignValW('Md3' ,DSQRT(PAR(9)))    

      CALL assignValW('Mha',AMASS(1))
      CALL assignValW('Mhb',AMASS(2))
      CALL assignValW('MHc',CMASS)

      CALL assignValW('muP',MUPSUSY)
      CALL assignValW('xvev',MUQ/LQ*DSQRT(ZS)) 

      DO i=1,3
       WRITE(name,fmt='(A2,I1)') 'Mh',i
       CALL assignValW(name,SMASS(i))
      ENDDO

      DO i=1,5
       WRITE(name,fmt='(A3,I1)') 'MNE',i
       CALL assignValW(name,MNEU(i))
       DO j=1,5
         WRITE(name,fmt='(A2,I1,I1)') 'Zn',i,j
         CALL assignValW(name,NEU(i,NORD(j)))
        ENDDO
      ENDDO

      CALL assignValW("Zh11",SCOMP(1,2))
      CALL assignValW("Zh12",SCOMP(1,1))
      CALL assignValW("Zh13",SCOMP(1,3))
      CALL assignValW("Zh21",SCOMP(2,2))
      CALL assignValW("Zh22",SCOMP(2,1))
      CALL assignValW("Zh23",SCOMP(2,3))
      CALL assignValW("Zh31",SCOMP(3,2))
      CALL assignValW("Zh32",SCOMP(3,1))
      CALL assignValW("Zh33",SCOMP(3,3))

      CALL assignValW('Za11',SINB*PCOMP(1,1))
      CALL assignValW('Za12',COSB*PCOMP(1,1))
      CALL assignValW('Za13',PCOMP(1,2))
      CALL assignValW('Za21',SINB*PCOMP(2,1))
      CALL assignValW('Za22',COSB*PCOMP(2,1))
      CALL assignValW('Za23',PCOMP(2,2)) 

      CALL assignValW('MSl1',MSL1)
      CALL assignValW('MSl2',MSL2)
      CALL assignValW('Zl11',CSL)
      CALL assignValW('Zl12',SSL)
      CALL assignValW('Zl21',-SSL)
      CALL assignValW('Zl22',CSL)

      CALL assignValW('MSb1',MSB1)
      CALL assignValW('MSb2',MSB2)
      CALL assignValW('Zb11',CSB)
      CALL assignValW('Zb12',SSB)
      CALL assignValW('Zb21',-SSB)
      CALL assignValW('Zb22',CSB)

      CALL assignValW('MSt1',MST1)
      CALL assignValW('MSt2',MST2)
      CALL assignValW('Zt11',CST)
      CALL assignValW('Zt12',SST)
      CALL assignValW('Zt21',-SST)
      CALL assignValW('Zt22',CST)


      CALL assignValW('Zu11',UU(1,1))
      CALL assignValW('Zu12',UU(1,2))
      CALL assignValW('Zu21',UU(2,1))
      CALL assignValW('Zu22',UU(2,2))

      CALL assignValW('Zv11',VV(1,1))
      CALL assignValW('Zv12',VV(1,2))
      CALL assignValW('Zv21',VV(2,1))
      CALL assignValW('Zv22',VV(2,2))

      CALL assignValW('MSeL',MLL)
      CALL assignValW('MSeR',MLR)
      CALL assignValW('MSmL',MLL)
      CALL assignValW('MSmR',MLR)
      CALL assignValW('MSne',MNL)
      CALL assignValW('MSnm',MNL)
      CALL assignValW('MSnl',MSNT)
      CALL assignValW('MSuL',MUL)
      CALL assignValW('MSuR',MUR)
      CALL assignValW('MSdL',MDL)
      CALL assignValW('MSdR',MDR)
      CALL assignValW('MScL',MUL)
      CALL assignValW('MScR',MUR)
      CALL assignValW('MSsL',MDL)
      CALL assignValW('MSsR',MDR)
      CALL assignValW('MSG',MGL)
      CALL assignValW('MC1',MCHA(1))
      CALL assignValW('MC2',MCHA(2))

*   Improved Higgs potential

       CALL assignValW("la1",PL(1))
       CALL assignValW("la2",PL(2))
       CALL assignValW("la3",PL(3))
       CALL assignValW("la4",PL(4))
       CALL assignValW("la5",PL(5))
       CALL assignValW("la6",PL(6))
       CALL assignValW("la7",PL(7))
       CALL assignValW("la1s",PK(1))
       CALL assignValW("la2s",PK(2))
       CALL assignValW("la3s",PK(3))
       CALL assignValW("la4s",PK(4)) 
       CALL assignValW("la5s",PK(5))
       CALL assignValW("la6s",PK(6))
       CALL assignValW("la7s",PK(7))
       CALL assignValW("la8s",PK(8))
       CALL assignValW("aa1",PA(1))
       CALL assignValW("aa2",PA(2))
       CALL assignValW("aa3",PA(3))
       CALL assignValW("aa4",PA(4))
       CALL assignValW("aa5",PA(5)) 
       CALL assignValW("aa6",PA(6))
       CALL assignValW("B1",PB(1))
       CALL assignValW("B2",PB(2))
       CALL assignValW("X",PX)
       CALL assignValW("dMb",DELMB)
 
*   Sorting sparticles

      err=sortOddParticles(mess)
      IF(mess.ne.'~o1'.AND.(M32.GT.DABS(MNEU(1)).OR.GRFLAG.EQ.0)) THEN
        OMG=-1d0
        PROB(30)=-1d0
        RETURN
      ENDIF
      IF(err.ne.0) THEN
        OMG=-2d0
        PROB(30)=-2d0
        RETURN
      ENDIF

      call ModelConstIni(2,WIDTH,err)

*   Computing relic density

      call forceUG(1)    ! set to 1/0 for unitary gauge on/off

      IF(MOFLAG.EQ.0 .OR. MOFLAG.EQ.2 .OR. MOFLAG.EQ.4 .OR. MOFLAG.EQ.6)
     . THEN
       fast=1
      ELSE
       fast=0
      ENDIF
      IF(MOFLAG.EQ.0 .OR. MOFLAG.EQ.1 .OR. MOFLAG.EQ.4 .OR. MOFLAG.EQ.5)
     . THEN
       Beps=1d-3
      ELSE
       Beps=1d-6
      ENDIF
      IF(MOFLAG.EQ.0 .OR. MOFLAG.EQ.1 .OR. MOFLAG.EQ.2 .OR. MOFLAG.EQ.3)
     . THEN
       WW=0
      ELSE
       WW=1
      ENDIF
      call setVVdecay(WW,0) ! set to 1/0 for vitual decays on/off
      OMG=darkOmega(Xf,fast,Beps,err)

      IF(M32.LT.DABS(MNEU(1)) .AND. GRFLAG.EQ.1)THEN
       OMG=OMG*M32/DABS(MNEU(1))
      ENDIF

      IF(OMG.GT.0d0)THEN
       IF(OMGFLAG.GT.0)THEN
        PROB(30)=DDIM(OMG/OMGMAX,1d0)-DDIM(1d0,OMG/OMGMIN)
       ELSE
        PROB(30)=DDIM(OMG/OMGMAX,1d0)
       ENDIF
      ELSE
       PROB(30)=-3d0
       RETURN
      ENDIF

*  (in)Direct detection constraints

      csPsi=0d0
      csNsi=0d0
      csPsd=0d0
      csNsd=0d0
      vcsll=0d0
      vcsbb=0d0
      sigmaV=0d0
      DO I=1,NBIN
       dNdx(I)=0d0
      ENDDO
      IF(M32.LT.DABS(MNEU(1)) .AND. GRFLAG.EQ.1) RETURN
      IF (IABS(OMGFLAG).EQ.1) RETURN
      IF (IABS(OMGFLAG).EQ.3) GOTO 1

*  Computing WIMP-Nucleon cross sections
*  Muq/Mdq=0.553d0, Msq/Mdq=18.9d0

      CALL calcScalarQuarkFF(0.553d0,18.9d0,sigmaPiN,sigmaS)
      err=nucleonAmplitudes(CDM1,pA0,pA5,nA0,nA5)
      Nmass=0.939d0
      SCcoeff=4d0/PI*3.8937966d8*(Nmass*Mcdm/(Nmass+Mcdm))**2
      csPsi=SCcoeff*pA0(1)**2
      csNsi=SCcoeff*nA0(1)**2
      csPsd=3*SCcoeff*pA5(1)**2
      csNsd=3*SCcoeff*nA5(1)**2
      IF( pA0(1)*nA0(1) .lt. 0d0) csNsi=-csNsi
      IF( pA5(1)*nA5(1) .lt. 0d0) csNsd=-csNsd
* New June 2019
      IF(OMGFLAG.LT.0 .AND.OMG.LT.OMGMAX) THEN
       scal=OMG/OMGMAX
      ELSE
       scal=1d0
      ENDIF
* End New

      lim=MIN(PandaX_SI(DABS(MNEU(1))),LUX_SI(DABS(MNEU(1))),
     .        XENON_SI(DABS(MNEU(1))),CRESST_SI(DABS(MNEU(1))),
     .        DarkSide50_SI(DABS(MNEU(1))))
      PROB(31)=DDIM(scal*DABS(csPsi)/lim,1d0)

      lim=MIN(LUX_SDn(DABS(MNEU(1))),XENON_SDn(DABS(MNEU(1))))
      PROB(61)=DDIM(scal*DABS(csNsd)/lim,1d0)

      lim=MIN(LUX_SDp(DABS(MNEU(1))),XENON_SDp(DABS(MNEU(1))),
     .        PICO60_SDp(DABS(MNEU(1))))
      PROB(62)=DDIM(scal*DABS(csPsd)/lim,1d0)

      IF (IABS(OMGFLAG).EQ.2) RETURN

*  Computing indirect detection rate
 1    sigmaV=calcSpectrum(0,tab,NULL,NULL,NULL,NULL,NULL,err)
      IF (err.NE.0) sigmaV=0d0
      IF (sigmaV.NE.0d0) THEN
       DO I=1,NBIN
        dNdx(I)=zInterp(DLOG(10d0)*x(I),tab)*DLOG(10d0)
       ENDDO
      ENDIF
      I=1
 123  ok=vSigmaCh(I,weight,pdg,process)
      if(ok.eq.1) then
      if((abs(pdg(3)).eq.15).and.(abs(pdg(4))).eq.15)vcsll=sigmaV*weight
      if((abs(pdg(3)).eq.5).and.(abs(pdg(4))).eq.5)vcsbb=sigmaV*weight
      I=I+1
      goto 123 
      endif

      END


      DOUBLE PRECISION FUNCTION XENON_SI(M)

*  Limit on the spin independent WIMP-nucleon cross-section,
*  from XENON_1T [arXiv:1805.12562]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=89)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/6.0000d0,6.2104d0,6.3384d0,6.4593d0,6.5862d0,6.7472d0,
     .6.9091d0,7.0580d0,7.2116d0,7.3882d0,7.5727d0,7.7317d0,7.9401d0,
     .8.1289d0,8.3779d0,8.6151d0,8.8929d0,9.1846d0,9.4770d0,9.8199d0,
     .1.0176d1,1.0532d1,1.0912d1,1.1394d1,1.1903d1,1.2428d1,1.2951d1,
     .1.3611d1,1.4311d1,1.5040d1,1.5887d1,1.6851d1,1.7841d1,1.9021d1,
     .2.0272d1,2.1646d1,2.3109d1,2.4731d1,2.6635d1,2.8594d1,3.0698d1,
     .3.3191d1,3.5887d1,3.8527d1,4.1715d1,4.4963d1,4.8381d1,5.1910d1,
     .5.5771d1,5.9829d1,6.4003d1,6.8913d1,7.3745d1,7.9476d1,8.5322d1,
     .9.2004d1,9.9101d1,1.0616d2,1.1415d2,1.2255d2,1.3167d2,1.4111d2,
     .1.5120d2,1.6169d2,1.7353d2,1.8629d2,2.0038d2,2.1508d2,2.3120d2,
     .2.4747d2,2.6567d2,2.8418d2,3.0622d2,3.2873d2,3.5542d2,3.8431d2,
     .4.1197d2,4.4208d2,4.7460d2,5.0769d2,5.4512d2,5.8435d2,6.2717d2,
     .6.7811d2,7.3048d2,7.8189d2,8.4136d2,9.0079d2,1.0000d3/

      DATA Y/2.5000d-8,1.8669d-8,1.5484d-8,1.2790d-8,1.0624d-8,
     .8.5913d-9,7.0589d-9,5.8275d-9,4.8561d-9,4.0277d-9,3.3407d-9,
     .2.7708d-9,2.3259d-9,1.9421d-9,1.6108d-9,1.3361d-9,1.0977d-9,
     .9.0209d-10,7.5513d-10,6.2708d-10,5.1979d-10,4.3842d-10,3.7105d-10,
     .3.0775d-10,2.5528d-10,2.1599d-10,1.8550d-10,1.5676d-10,1.3268d-10,
     .1.1422d-10,9.7939d-11,8.4727d-11,7.4249d-11,6.5352d-11,5.8370d-11,
     .5.3340d-11,4.8918d-11,4.6191d-11,4.4386d-11,4.3173d-11,4.2646d-11,
     .4.2760d-11,4.3510d-11,4.4585d-11,4.6475d-11,4.8702d-11,5.1189d-11,
     .5.3904d-11,5.6923d-11,5.9907d-11,6.3563d-11,6.7311d-11,7.1611d-11,
     .7.6081d-11,8.0916d-11,8.6440d-11,9.2432d-11,9.7879d-11,1.0500d-10,
     .1.1208d-10,1.1947d-10,1.2682d-10,1.3546d-10,1.4424d-10,1.5468d-10,
     .1.6459d-10,1.7657d-10,1.8908d-10,2.0274d-10,2.1713d-10,2.3230d-10,
     .2.4792d-10,2.6716d-10,2.8624d-10,3.0949d-10,3.3371d-10,3.5560d-10,
     .3.8270d-10,4.1005d-10,4.4080d-10,4.7298d-10,5.0753d-10,5.4244d-10,
     .5.9060d-10,6.3378d-10,6.8016d-10,7.3066d-10,7.8997d-10,8.8000d-10/

      XENON_SI=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        XENON_SI=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION XENON_SDn(M)

*  Limit on the spin dependent WIMP-neutron cross-section,
*  from XENON_1T [arXiv:1902.03234]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=40)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/6.0000d0,6.4028d0,6.8456d0,7.2649d0,7.7434d0,8.2693d0,
     .8.9346d0,9.6939d0,1.0645d1,1.1770d1,1.3177d1,1.4903d1,1.7024d1,
     .1.9608d1,2.2670d1,2.6411d1,3.0770d1,3.5848d1,4.1765d1,4.8657d1,
     .5.6905d1,6.6044d1,7.6650d1,8.9985d1,1.0484d2,1.2167d2,1.4175d2,
     .1.6515d2,1.9240d2,2.2330d2,2.5917d2,3.0194d2,3.5043d2,4.1140d2,
     .4.7746d2,5.5415d2,6.5055d2,7.5503d2,8.7628d2,1.0000d3/
      DATA Y/3.5260d-3,2.0887d-3,1.2241d-3,7.2632d-4,4.3097d-4,
     .2.5572d-4,1.5401d-4,9.4940d-5,5.8975d-5,3.7663d-5,2.4547d-5,
     .1.6644d-5,1.1653d-5,8.9305d-6,7.3100d-6,6.5299d-6,6.3537d-6,
     .6.7358d-6,7.5123d-6,8.4175d-6,9.5372d-6,1.0824d-5,1.2298d-5,
     .1.4145d-5,1.6117d-5,1.8546d-5,2.1334d-5,2.4592d-5,2.8317d-5,
     .3.2592d-5,3.7548d-5,4.3236d-5,5.0083d-5,5.8568d-5,6.7903d-5,
     .7.8400d-5,9.2140d-5,1.0664d-4,1.2429d-4,1.4500d-4/

      XENON_SDn=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        XENON_SDn=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION XENON_SDp(M)

*  Limit on the spin dependent WIMP-proton cross-section,
*  from XENON_1T [arXiv:1902.03234]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=39)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/6.0000d0,6.5361d0,6.9402d0,7.3797d0,7.9220d0,8.5613d0,
     .9.2990d0,1.0131d1,1.1254d1,1.2510d1,1.4073d1,1.6026d1,1.8390d1,
     .2.1427d1,2.4870d1,2.8978d1,3.3764d1,3.9340d1,4.5838d1,5.3613d1,
     .6.2467d1,7.3623d1,8.6111d1,1.0110d2,1.1780d2,1.3831d2,1.6239d2,
     .1.9066d2,2.2215d2,2.5884d2,3.0159d2,3.5006d2,4.0632d2,4.7706d2,
     .5.5373d2,6.4272d2,7.5462d2,8.7589d2,1.0000d3/
      DATA Y/1.2000d-1,5.3518d-2,3.0855d-2,1.8298d-2,1.0699d-2,
     .6.2561d-3,3.7629d-3,2.2633d-3,1.4265d-3,9.3240d-4,6.1732d-4,
     .4.3045d-4,3.1536d-4,2.4308d-4,2.0810d-4,1.9667d-4,2.0231d-4,
     .2.2650d-4,2.5359d-4,2.8428d-4,3.2022d-4,3.6610d-4,4.0988d-4,
     .4.7205d-4,5.3240d-4,6.0867d-4,7.0098d-4,8.0730d-4,9.2975d-4,
     .1.0708d-3,1.2332d-3,1.4202d-3,1.6426d-3,1.9324d-3,2.2315d-3,
     .2.5700d-3,3.0396d-3,3.5064d-3,4.0000d-3/

      XENON_SDp=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        XENON_SDp=1d36*(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION PandaX_SI(M)

*  Limit on the spin independent WIMP-nucleon cross-section,
*  from PandaX [arXiv:1607.07400]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=92)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/4.9961d0,5.0730d0,5.1498d0,5.2302d0,5.3197d0,5.4111d0,
     .5.5000d0,5.5982d0,5.6952d0,5.7899d0,5.8948d0,6.0092d0,6.1416d0,
     .6.2740d0,6.4202d0,6.5615d0,6.7156d0,6.8700d0,7.0265d0,7.2181d0,
     .7.4017d0,7.6059d0,7.8082d0,8.0320d0,8.2647d0,8.5166d0,8.7967d0,
     .9.1045d0,9.4388d0,9.7992d0,1.0159d1,1.0551d1,1.0978d1,1.1383d1,
     .1.1927d1,1.2500d1,1.3095d1,1.3730d1,1.4436d1,1.5195d1,1.6074d1,
     .1.7050d1,1.8106d1,1.9317d1,2.0606d1,2.2100d1,2.3737d1,2.5590d1,
     .2.7843d1,3.0202d1,3.3021d1,3.6104d1,3.9306d1,4.2792d1,4.6588d1,
     .5.0936d1,5.5606d1,6.0630d1,6.6008d1,7.2169d1,7.8570d1,8.5177d1,
     .9.2642d1,1.0053d2,1.0898d2,1.1815d2,1.2808d2,1.3932d2,1.5158d2,
     .1.6457d2,1.7841d2,1.9387d2,2.1040d2,2.2746d2,2.4642d2,2.6811d2,
     .2.9003d2,3.1395d2,3.4179d2,3.6987d2,3.9999d2,4.3328d2,4.6978d2,
     .5.0875d2,5.5245d2,5.9637d2,6.4598d2,7.0083d2,7.5980d2,8.2020d2,
     .8.8164d2,9.5360d2/
      DATA Y/9.2058d-7,7.6741d-7,6.3394d-7,5.2846d-7,4.4053d-7,
     .3.6391d-7,3.0336d-7,2.5289d-7,2.1081d-7,1.7574d-7,1.4705d-7,
     .1.2212d-7,1.0088d-7,8.4096d-8,6.9469d-8,5.7910d-8,4.8275d-8,
     .4.0026d-8,3.3547d-8,2.7712d-8,2.3101d-8,1.9126d-8,1.6054d-8,
     .1.3383d-8,1.1258d-8,9.4705d-9,7.9670d-9,6.6414d-9,5.5000d-9,
     .4.5000d-9,3.8000d-9,3.2000d-9,2.6500d-9,2.2500d-9,1.9000d-9,
     .1.5943d-9,1.3526d-9,1.1516d-9,9.8404d-10,8.4397d-10,7.2512d-10,
     .6.2584d-10,5.4685d-10,4.7814d-10,4.2235d-10,3.7635d-10,3.4003d-10,
     .3.1014d-10,2.8662d-10,2.7154d-10,2.5938d-10,2.5149d-10,2.4804d-10,
     .2.4892d-10,2.5006d-10,2.5529d-10,2.6385d-10,2.7223d-10,2.8520d-10,
     .2.9946d-10,3.1458d-10,3.3311d-10,3.5290d-10,3.7501d-10,3.9869d-10,
     .4.2427d-10,4.5096d-10,4.8082d-10,5.1261d-10,5.4563d-10,5.8115d-10,
     .6.2070d-10,6.6702d-10,7.1843d-10,7.7541d-10,8.3732d-10,9.0173d-10,
     .9.7180d-10,1.0499d-9,1.1314d-9,1.2190d-9,1.3196d-9,1.4321d-9,
     .1.5553d-9,1.7010d-9,1.8575d-9,2.0234d-9,2.2148d-9,2.4239d-9,
     .2.6522d-9,2.9208d-9,3.2220d-9/

      PandaX_SI=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        PandaX_SI=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LUX_SI(M)

*  Limit on the spin independent WIMP-nucleon cross-section,
*  from LUX [arXiv:1608.07648]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=99)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/5.4616d0,5.5230d0,5.5883d0,5.6779d0,5.7729d0,5.8659d0,
     .5.9617d0,6.0578d0,6.1569d0,6.2568d0,6.3674d0,6.5140d0,6.6628d0,
     .6.8195d0,6.9737d0,7.1240d0,7.3222d0,7.5346d0,7.7384d0,7.9899d0,
     .8.3084d0,8.6400d0,9.0234d0,9.4196d0,9.8401d0,1.0260d1,1.0731d1,
     .1.1224d1,1.1803d1,1.2410d1,1.3206d1,1.4109d1,1.5282d1,1.6640d1,
     .1.8303d1,2.0054d1,2.2301d1,2.5339d1,2.8659d1,3.2717d1,3.7766d1,
     .4.3862d1,5.1243d1,6.0149d1,7.0327d1,8.2227d1,9.6141d1,1.1241d2,
     .1.2983d2,1.4894d2,1.7144d2,1.9697d2,2.2716d2,2.5944d2,2.9720d2,
     .3.3706d2,3.8425d2,4.3932d2,5.0040d2,5.6947d2,6.5254d2,7.4578d2,
     .8.4754d2,9.6815d2,1.1060d3,1.2631d3,1.4443d3,1.6572d3,1.8944d3,
     .2.1489d3,2.4514d3,2.8016d3,3.1990d3,3.6560d3,4.1763d3,4.7359d3,
     .5.3979d3,6.1298d3,7.0323d3,7.9989d3,9.0461d3,1.0334d4,1.1807d4,
     .1.3414d4,1.5301d4,1.7568d4,1.9967d4,2.2583d4,2.5793d4,2.9468d4,
     .3.3407d4,3.8028d4,4.3233d4,4.9370d4,5.6294d4,6.4376d4,7.3073d4,
     .8.3439d4,9.5764d4/
      DATA Y/8.7557d-7,7.0640d-7,5.6991d-7,4.5245d-7,3.6307d-7,
     .2.9292d-7,2.3633d-7,1.9066d-7,1.5300d-7,1.2278d-7,9.9055d-8,
     .7.9063d-8,6.3787d-8,5.1462d-8,4.1076d-8,3.3139d-8,2.6736d-8,
     .2.1570d-8,1.7403d-8,1.4040d-8,1.1327d-8,9.1387d-9,7.2942d-9,
     .5.8849d-9,4.7478d-9,3.8305d-9,3.1092d-9,2.5202d-9,2.0332d-9,
     .1.6404d-9,1.3377d-9,1.0830d-9,8.9921d-10,7.4657d-10,6.1715d-10,
     .5.1648d-10,4.3820d-10,3.7100d-10,3.2163d-10,2.8051d-10,2.5193d-10,
     .2.2884d-10,2.2296d-10,2.2527d-10,2.3329d-10,2.4766d-10,2.6208d-10,
     .2.8393d-10,3.0905d-10,3.4542d-10,3.8158d-10,4.2559d-10,4.8292d-10,
     .5.4193d-10,6.0874d-10,6.9144d-10,7.8903d-10,8.9594d-10,1.0213d-9,
     .1.1812d-9,1.3503d-9,1.5340d-9,1.7598d-9,2.0227d-9,2.2995d-9,
     .2.6163d-9,3.0013d-9,3.4343d-9,3.9009d-9,4.4299d-9,5.0507d-9,
     .5.7378d-9,6.5244d-9,7.4947d-9,8.5192d-9,9.6949d-9,1.1157d-8,
     .1.2656d-8,1.4477d-8,1.6406d-8,1.8901d-8,2.1493d-8,2.4423d-8,
     .2.7992d-8,3.2189d-8,3.6818d-8,4.1737d-8,4.8049d-8,5.4667d-8,
     .6.2130d-8,7.0651d-8,8.1338d-8,9.2339d-8,1.0627d-7,1.2106d-7,
     .1.3908d-7,1.5797d-7,1.8127d-7,2.0872d-7/

      LUX_SI=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        LUX_SI=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LUX_SDn(M)

*  Limit on the spin dependent WIMP-neutron cross-section,
*  from LUX [arXiv:1602.03489]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=65)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/3.5089d0,3.5559d0,3.6117d0,3.6944d0,3.8005d0,3.9322d0,
     .4.0587d0,4.2124d0,4.3938d0,4.5995d0,4.8633d0,5.1178d0,5.3854d0,
     .5.7481d0,6.1375d0,6.6137d0,7.2497d0,7.9793d0,8.9073d0,1.0198d1,
     .1.1803d1,1.3756d1,1.6562d1,2.0270d1,2.4841d1,3.0667d1,3.8134d1,
     .4.7077d1,5.8117d1,7.1746d1,8.8571d1,1.0934d2,1.3401d2,1.6424d2,
     .2.0128d2,2.4669d2,3.0454d2,3.7323d2,4.5743d2,5.6122d2,6.8584d2,
     .8.3819d2,1.0245d3,1.2556d3,1.5501d3,1.8998d3,2.3515d3,2.8743d3,
     .3.5484d3,4.3488d3,5.3686d3,6.6276d3,8.1819d3,1.0101d4,1.2469d4,
     .1.5394d4,1.8866d4,2.3122d4,2.8337d4,3.4729d4,4.2563d4,5.2545d4,
     .6.4398d4,7.8924d4,9.7433d4/
      DATA Y/5.8308d1,3.0609d1,1.6417d1,8.8053d0,4.6223d0,2.4792d0,
     .1.3015d0,6.9803d-1,3.7439d-1,2.0080d-1,1.0576d-1,5.7766d-2,
     .3.0983d-2,1.6837d-2,9.3041d-3,4.9902d-3,2.7610d-3,1.5644d-3,
     .9.1432d-4,5.5074d-4,3.3725d-4,2.1965d-4,1.5809d-4,1.2213d-4,
     .1.0363d-4,9.5468d-5,9.3017d-5,1.0144d-4,1.2024d-4,1.2829d-4,
     .1.4925d-4,1.7713d-4,2.1255d-4,2.5532d-4,3.0719d-4,3.6834d-4,
     .4.4604d-4,5.4319d-4,6.6426d-4,8.1337d-4,9.8598d-4,1.1972d-3,
     .1.4534d-3,1.7725d-3,2.1856d-3,2.6573d-3,3.2858d-3,3.9900d-3,
     .4.9085d-3,6.0098d-3,7.4498d-3,9.2291d-3,1.1438d-2,1.4175d-2,
     .1.7569d-2,2.1777d-2,2.6812d-2,3.2993d-2,4.0385d-2,4.9253d-2,
     .5.9800d-2,7.4685d-2,9.1296d-2,1.1143d-1,1.3815d-1/

      LUX_SDn=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        LUX_SDn=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION LUX_SDp(M)

*  Limit on the spin dependent WIMP-proton cross-section,
*  from LUX [arXiv:1602.03489]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=63)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/3.7146d0,3.8006d0,3.9158d0,4.0283d0,4.1475d0,4.3036d0,
     .4.5080d0,4.7211d0,4.9366d0,5.2311d0,5.5566d0,5.9213d0,6.3241d0,
     .6.8525d0,7.5330d0,8.3292d0,9.4195d0,1.0801d1,1.2426d1,1.4850d1,
     .1.7937d1,2.1940d1,2.6886d1,3.3188d1,4.1115d1,5.0936d1,6.2886d1,
     .7.7611d1,9.5801d1,1.1740d2,1.4387d2,1.7630d2,2.1605d2,2.6475d2,
     .3.2444d2,3.9759d2,4.8722d2,5.9706d2,7.3167d2,8.9925d2,1.0988d3,
     .1.3465d3,1.6500d3,2.0220d3,2.4959d3,3.0587d3,3.7482d3,4.5932d3,
     .5.6288d3,6.9480d3,8.5144d3,1.0434d4,1.2879d4,1.5783d4,1.9341d4,
     .2.3702d4,2.8927d4,3.5336d4,4.3302d4,5.3452d4,6.5641d4,8.0269d4,
     .9.8022d4/
      DATA Y/4.2591d2,2.2492d2,1.2054d2,6.4599d1,3.4620d1,1.8553d1,
     .9.7316d0,5.2154d0,2.7950d0,1.4979d0,8.2021d-1,4.4912d-1,2.4593d-1,
     .1.3759d-1,7.6979d-2,4.3762d-2,2.5462d-2,1.5474d-2,9.7035d-3,
     .6.6681d-3,4.8470d-3,3.7559d-3,3.0776d-3,2.8994d-3,3.0301d-3,
     .3.0744d-3,3.3700d-3,3.8407d-3,4.3656d-3,5.1836d-3,6.2232d-3,
     .7.5406d-3,9.0686d-3,1.0909d-2,1.3217d-2,1.6280d-2,2.0104d-2,
     .2.4668d-2,3.0377d-2,3.7481d-2,4.5532d-2,5.5336d-2,6.7269d-2,
     .8.1867d-2,9.9782d-2,1.2150d-1,1.4781d-1,1.8025d-1,2.1910d-1,
     .2.6932d-1,3.2924d-1,4.0028d-1,4.9247d-1,6.0128d-1,7.3156d-1,
     .8.9410d-1,1.0974d0,1.3330d0,1.6186d0,2.0243d0,2.4853d0,3.0205d0,
     .3.7224d0/

      LUX_SDp=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        LUX_SDp=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION PICO60_SDp(M)

*  Limit on the spin dependent WIMP-proton cross-section,
*  from PICO-60 [arXiv:1902.04031]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=64)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/2.9000d0,2.9323d0,2.9623d0,2.9924d0,3.0182d0,3.0464d0,
     .3.1114d0,3.1734d0,3.2409d0,3.3053d0,3.3764d0,3.4443d0,3.5187d0,
     .3.5896d0,3.6607d0,3.7719d0,3.8998d0,4.0398d0,4.2004d0,4.3686d0,
     .4.5767d0,4.7915d0,5.0232d0,5.3188d0,5.7249d0,6.1316d0,6.6771d0,
     .7.2787d0,8.1158d0,9.0561d0,1.0247d1,1.1589d1,1.3293d1,1.5244d1,
     .1.7481d1,2.0046d1,2.2989d1,2.6363d1,3.0440d1,3.4907d1,4.0582d1,
     .4.6538d1,5.3369d1,6.2045d1,7.1152d1,8.2719d1,9.4860d1,1.0878d2,
     .1.2390d2,1.4111d2,1.6405d2,1.8813d2,2.1574d2,2.4741d2,2.8567d2,
     .3.2759d2,3.7825d2,4.3081d2,4.9067d2,5.7044d2,6.4526d2,7.5017d2,
     .8.5440d2,1.0000d3/
      DATA Y/3.0000d-1,2.0000d-1,1.4471d-1,1.0310d-1,7.6774d-2,
     .5.6631d-2,4.1337d-2,3.0332d-2,2.1610d-2,1.5857d-2,1.1297d-2,
     .8.2894d-3,5.9058d-3,4.3335d-3,3.1798d-3,2.3679d-3,1.7375d-3,
     .1.2379d-3,9.2181d-4,6.7639d-4,5.1117d-4,3.8630d-4,2.9194d-4,
     .2.1743d-4,1.6429d-4,1.2600d-4,9.7747d-5,7.7484d-5,6.1129d-5,
     .4.9674d-5,4.1692d-5,3.5771d-5,3.1976d-5,2.9369d-5,2.7594d-5,
     .2.6680d-5,2.6096d-5,2.6062d-5,2.6420d-5,2.7161d-5,2.8081d-5,
     .2.9657d-5,3.1432d-5,3.3548d-5,3.6514d-5,3.9671d-5,4.3423d-5,
     .4.7640d-5,5.2495d-5,5.8264d-5,6.5512d-5,7.3859d-5,8.3063d-5,
     .9.3318d-5,1.0456d-4,1.1763d-4,1.3366d-4,1.5206d-4,1.7410d-4,
     .1.9682d-4,2.2614d-4,2.5581d-4,2.9247d-4,3.3335d-4/

      PICO60_SDp=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        PICO60_SDp=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION CRESST_SI(M)

*  Limit on the spin independent WIMP-nucleon cross-section,
*  from CRESST-II [arXiv:1509.01515]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=60)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/5.0000d-1,5.3990d-1,5.6315d-1,5.8848d-1,6.1495d-1,
     .6.4261d-1,6.7399d-1,7.0689d-1,7.4413d-1,7.8332d-1,8.2761d-1,
     .8.8084d-1,9.4439d-1,1.0088d0,1.0816d0,1.1597d0,1.2479d0,1.3379d0,
     .1.4345d0,1.5493d0,1.6733d0,1.8006d0,1.9519d0,2.1081d0,2.2685d0,
     .2.4322d0,2.6172d0,2.8061d0,3.0085d0,3.2732d0,3.5030d0,3.7626d0,
     .4.0489d0,4.5197d0,4.8636d0,5.2337d0,5.6526d0,6.0827d0,6.5575d0,
     .7.0694d0,7.6073d0,8.1862d0,8.8090d0,9.4793d0,1.0201d1,1.0977d1,
     .1.1812d1,1.2711d1,1.3678d1,1.4773d1,1.5897d1,1.7106d1,1.8408d1,
     .1.9954d1,2.1473d1,2.3106d1,2.4864d1,2.6805d1,2.8582d1,3.0000d1/
      DATA Y/1.8393d0,1.0645d0,6.7862d-1,4.2193d-1,2.6233d-1,1.6310d-1,
     .1.0141d-1,6.6283d-2,4.3325d-2,2.9771d-2,2.0976d-2,1.5153d-2,
     .1.1224d-2,8.5244d-3,6.6380d-3,5.3000d-3,4.2317d-3,3.5520d-3,
     .2.9815d-3,2.5660d-3,2.3217d-3,2.0488d-3,1.8079d-3,1.5560d-3,
     .1.2738d-3,9.4352d-4,6.4835d-4,4.4553d-4,3.0615d-4,1.8334d-4,
     .1.4099d-4,1.2134d-4,9.9338d-5,8.1323d-5,7.3580d-5,6.6575d-5,
     .5.8749d-5,5.1842d-5,4.4618d-5,4.1392d-5,3.8400d-5,3.7452d-5,
     .3.5624d-5,3.4744d-5,3.3886d-5,3.3886d-5,3.4744d-5,3.4744d-5,
     .3.6526d-5,3.8400d-5,3.9373d-5,4.1392d-5,4.2441d-5,4.5748d-5,
     .4.6907d-5,4.9313d-5,5.0562d-5,4.5748d-5,3.5624d-5,2.9560d-5/

      CRESST_SI=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        CRESST_SI=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION DarkSide50_SI(M)

*  Limit on the spin independent WIMP-nucleon cross-section,
*  from DarkSide-50 [arXiv:1802.06994]
*  X = WIMP mass in GeV, Y = cross section in pb

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=14)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/1.8000d0,2.1317d0,2.5033d0,2.9464d0,3.4058d0,3.9468d0,
     .4.5983d0,5.3724d0,6.3551d0,7.5255d0,8.8755d0,1.0510d1,1.2395d1,
     .1.5000d1/
      DATA Y/1.3040d-5,9.0862d-6,6.6089d-6,5.3928d-6,2.9517d-6,
     .1.6005d-6,9.6411d-7,6.2724d-7,4.5366d-7,3.5509d-7,2.9509d-7,
     .2.7191d-7,2.5886d-7,2.5886d-7/

      DarkSide50_SI=1d99
      IF(M.LT.X(1) .OR. M.GE.X(N))RETURN
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LT.X(I+1)))THEN
        DarkSide50_SI=Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I))
        RETURN
       ENDIF
      ENDDO

      END


      SUBROUTINE PrintRelDen(ch)

      IMPLICIT NONE

#include "../sources/micromegas/include/micromegas.fh"

      CHARACTER*50 txt 

      INTEGER NBIN,OMGFLAG,MAFLAG,MOFLAG
      INTEGER ch,I,pdg(5),err

      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX,Xf
      DOUBLE PRECISION sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION MGL,MCHA(2),UU(2,2),VV(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION weightCh,alph_,omgfo,v,nngg,nngz,ccoeff
      DOUBLE PRECISION vcsnngg,vcsnngz

      EXTERNAL vcsnngg,vcsnngz

      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/SUSYSPEC/MGL,MCHA,UU,VV,MNEU,NEU
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd

      IF(OMGFLAG.EQ.0)THEN

       WRITE(ch,900) 4,"# Omega h^2 not computed (OMGFLAG=0)"
      
      ELSE

       IF(OMG.GT.0d0)THEN
        WRITE(ch,899) "#"
        WRITE(ch,899) "BLOCK ABUNDANCE"
        WRITE(ch,920) 1, Mcdm/Xf,'T_f[GeV]'
        WRITE(ch,920) 4, OMG,"omega h^2"
        WRITE(ch,920) 3, vSigma(Mcdm/Xf,1d-3,1,alph_)*2.9979d-26,
     .  'vSigma'
        WRITE(ch,899)'# contibutions to vSigma in percents'
        I=1
1111    err=vSigmaTch(I,weightCh, pdg,txt)
        if(err.eq.1 .and. weightCh.gt.1.E-2) then
         write(ch,921)6,I,pdg(1),pdg(2),pdg(3),pdg(4),100*weightCh,txt
         I=I+1  
         goto 1111
        endif
        omgfo=PRINTCHANNELS(Xf,1d-3,1.D-4,1,0)
        WRITE(ch,899)'# contibutions to 1/Omega in percents'
        I=1
1112    err=omegaCh(I,weightCh, pdg,txt)
        if(err.eq.1 .and. weightCh.gt.1.E-2) then
         write(ch,921)7,I,pdg(1),pdg(2),pdg(3),pdg(4),100*weightCh,txt
         I=I+1  
         goto 1112
        endif
        WRITE(ch,899) "#"
        WRITE(ch,899) "BLOCK LSP"
        WRITE(ch,920) 0, Mcdm,"LSP mass"
        WRITE(ch,920) 1, dabs(NEU(1,1)), "bino"
        WRITE(ch,920) 2, dabs(NEU(1,2)), "wino" 
        WRITE(ch,920) 3, dabs(NEU(1,3)), "higgsino2"
        WRITE(ch,920) 4, dabs(NEU(1,4)), "higgsino1"
        WRITE(ch,920) 5, dabs(NEU(1,5)), "singlino"
        IF(IABS(OMGFLAG).EQ.2 .OR. IABS(OMGFLAG).EQ.4)THEN
         WRITE(ch,899) "#"
         WRITE(ch,899) 'BLOCK NDMCROSSSECT'
         WRITE(ch,920) 1,csPsi*1d-36,"# csPsi [cm^2]"
         WRITE(ch,920) 2,csNsi*1d-36,"# csNsi [cm^2]"
         WRITE(ch,920) 3,csPsd*1d-36,"# csPsd [cm^2]"
         WRITE(ch,920) 4,csNsd*1d-36,"# csNsd [cm^2]"
         WRITE(ch,925)"# Values used for sigma_piN,sigma_S",
     .   " (strange content of the proton)"
         WRITE(ch,924) "#",sigmapiN,"sigma_piN"
         WRITE(ch,924) "#",sigmaS,"sigma_S"
         IF(OMG.LT.OMGMIN)THEN
          WRITE(ch,925) "# Cross sections not rescaled with",
     .    " the ratio relic density/Planck observed value"
         ENDIF
        ENDIF
        IF(IABS(OMGFLAG).EQ.3 .OR. IABS(OMGFLAG).EQ.4)THEN
         WRITE(ch,899) "#"
         WRITE(ch,899) 'BLOCK ANNIHILATION'
         v=Vrot/299792d0*1.5957691d0
         nngg=vcsnngg(v)*2.9979d-26*0.9117d0**2
         nngz=vcsnngz(v)*2.9979d-26*0.9117d0
         ccoeff=sigmaV/(sigmaV+nngg+nngz)
         WRITE(ch,920) 0,sigmaV+nngg+nngz,'sigmaV [cm^3/s]'
         write(ch,923)1,nngg/(sigmaV+nngg+nngz),2,22,22,'~o1,~o1 -> A,A'
         write(ch,923)2,nngz/(sigmaV+nngg+nngz),2,22,23,'~o1,~o1 -> A,Z'
         I=3
 1113    err=vSigmach(I-2,weightCh, pdg,txt)
         if(err.eq.1 .and. weightCh.gt.1d-3) then
         if(.not.(
     .              ((pdg(3).eq.22).and.(pdg(4).eq.22))
     .          .or.((pdg(3).eq.22).and.(pdg(4).eq.23))
     .          .or.((pdg(3).eq.23).and.(pdg(4).eq.22))
     .           )) write(ch,923) I,weightCh*ccoeff,2,pdg(3),pdg(4),txt
          I=I+1
          goto 1113
         endif
        ENDIF

       ELSE
        WRITE(ch,900) 4,"# Problem in micrOMEGAs"
       ENDIF

      ENDIF

 899  FORMAT(A)
 900  FORMAT(1X,I5,3X,A)
 920  FORMAT(1x,I3,1x, 1PE12.4,3x,'# ',1x,A) 
 921  FORMAT(1x,I3,1x,I3,1x,I9,1x,I9,1x,I9,1x,I9,1x,1PE12.4,3x,'# ',A)
 923  FORMAT(1x,I3,1x,1PE12.4,1x,I2,1x,I9,1x,I9,3x,'# ',A)
 924  FORMAT(A,1P,E16.8,1X,A)
 925  FORMAT(A,A)

      end
