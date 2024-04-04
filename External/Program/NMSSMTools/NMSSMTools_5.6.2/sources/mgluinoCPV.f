       SUBROUTINE MGLUINO_CPV(PAR)

c       Gluino mass
c       The gluino mass MGL, including alpha_S corrections, is stored in the 
c       common GLUSPEC.

       IMPLICIT NONE

       DOUBLE PRECISION PAR(*),Pi
       DOUBLE PRECISION MGL2,ALSMT,ALSQ,S2T,S2B,DMG,DMQ,NMB0,NMB1,R
       DOUBLE PRECISION Q2
       DOUBLE PRECISION GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
       DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
       DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
       DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
       DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
       DOUBLE PRECISION MGL
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/RENSCALE/Q2
      COMMON/EWPAR/GF,MZ,MW,g1,g2,alSMZ,S2TW,ALEMMZ
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/GLUSPEC/MGL


      Pi=4d0*DATAN(1d0)
      MGL=PAR(22)
      MGL2=MGL**2
      ALSMT=ALSMZ/(1d0+23d0/(12d0*PI)*ALSMZ*DLOG(MT**2/MZ**2))
      ALSQ=ALSMT/(1d0+7d0*ALSMT/(4d0*pi)*DLOG(Q2/MT**2))
      S2T=2d0*UT(1,1,1)
     .       *(UT(1,2,1)*DDCOS(PhiM3)-UT(1,2,2)*DDSIN(PhiM3))
      S2B=2d0*UB(1,1,1)
     .       *(UB(1,2,1)*DDCOS(PhiM3)-UB(1,2,2)*DDSIN(PhiM3))
      R=0d0 ! 1=MSbar, 0=DRbar scheme

*  gluon/gluino correction, eq.(22) of BMPZ,
*  plus possible shift to MSbar

      DMG= 3d0*(2d0*NMB0(MGL2,MGL2,0d0,Q2)-NMB1(MGL2,MGL2,0d0,Q2)
     .     -R/2d0)

*  quark-squark loops, eq.(D.44) of BMPZ

      DMQ= 2d0*(NMB1(MGL2,0d0,MSU2(1),Q2) + NMB1(MGL2,0d0,MSU2(2),Q2)
     .       + NMB1(MGL2,0d0,MSD2(1),Q2) + NMB1(MGL2,0d0,MSD2(2),Q2))
     .       + NMB1(MGL2,MT**2,MST2(1),Q2) + NMB1(MGL2,MT**2,MST2(2),Q2)
     .       + NMB1(MGL2,MB**2,MSB2(1),Q2) + NMB1(MGL2,MB**2,MSB2(2),Q2)
     .       + MT/MGL*S2T*(NMB0(MGL2,MT**2,MST2(1),Q2)
     .       - NMB0(MGL2,MT**2,MST2(2),Q2))
     .       + MB/MGL*S2B*(NMB0(MGL2,MB**2,MSB2(1),Q2)
     .       - NMB0(MGL2,MB**2,MSB2(2),Q2))

      MGL=MGL/(1d0-ALSQ/(2d0*PI)*(DMG-DMQ/2d0))

c      print*,'MGL',MGL

       RETURN
       END
