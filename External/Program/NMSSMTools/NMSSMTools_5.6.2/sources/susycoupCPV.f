      SUBROUTINE SUSYCOUP_CPV()

C         Trilinear couplings involving the SUSY particles 
C      - Chargino/neutralino - SM fermion - sfermion couplings are computed
c        and stored in the commons CHSFfCOUP (charginos) and NEUSFfCOUP 
c        (neutralinos).
c      - Chargino/neutralino - Higgs couplings are also computed and stored in 
c        the common HINOCOUP.

      IMPLICIT NONE

      INTEGER I,J,M

      DOUBLE PRECISION gg1,gg2,Ytau,DELT(2,2),XHG(5,6)
      DOUBLE PRECISION G1Q,G2Q,GQ,ALSQ
      DOUBLE PRECISION ZHU,ZHD,ZS,vuq,vdq,TANBQ
      DOUBLE PRECISION tanb,cosb,sinb,vu,vd
      DOUBLE PRECISION mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      DOUBLE PRECISION l,k,Alcos1,Akcos2,muq,nuq
      DOUBLE PRECISION Ytq,Ybq,MTOPQ,MBOTQ
      DOUBLE PRECISION phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      DOUBLE PRECISION MST2(2),UT(2,2,2),MSB2(2),UB(2,2,2),MSL2(2),
     . UTAU(2,2,2),MSNT2
      DOUBLE PRECISION MSU2(2),MSD2(2),MSE2(2),MSNE2,MSMU2(2),
     . UMU(2,2,2)
      DOUBLE PRECISION MCH2(2),U(2,2,2),V(2,2,2)
      DOUBLE PRECISION MNEU(5),NEU(5,5,2)
      DOUBLE PRECISION MHC,XC(2,2),MH0(5),XH(5,5),MA2
      DOUBLE PRECISION COCHSTbL(2,2,2),COCHSTbR(2,2,2),
     . COCHSBtL(2,2,2),COCHSBtR(2,2,2),COCHSLnL(2,2,2),
     . COCHSNlL(2,2),COCHSNlR(2,2),COCHSUdL(2,2,2),COCHSUdR(2,2,2),
     . COCHSDuL(2,2,2),COCHSDuR(2,2,2),COCHSEnL(2,2,2),
     . COCHSNeL(2,2),COCHSNeR(2,2)
      DOUBLE PRECISION CONESTtL(5,2,2),CONESTtR(5,2,2),
     . CONESBbL(5,2,2),CONESBbR(5,2,2),CONESLlL(5,2,2),CONESLlR(5,2,2),
     . CONESNnL(5,2),CONESUuL(5,2,2),CONESUuR(5,2,2),
     . CONESDdL(5,2,2),CONESDdR(5,2,2),CONESEeL(5,2,2),
     . CONESEeR(5,2,2)
      DOUBLE PRECISION COH0CH(5,2,2,2),COH0NEU(5,5,5,2),
     . COHPNEUCHM(2,5,2,2),COHMNEUCHP(2,5,2,2)
      DOUBLE PRECISION DDCOS,DDSIN

      COMMON/QGAUGE/G1Q,G2Q,GQ,ALSQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,vuq,vdq,TANBQ
      COMMON/TBPAR/tanb,cosb,sinb,vu,vd
      COMMON/SMFERM/mt,mb,mtau,mmu,mel,MS,MC,MBP,MPI,MSTRANGE
      COMMON/QPAR/l,k,Alcos1,Akcos2,muq,NUQ
      COMMON/QQUARK/Ytq,Ybq,MTOPQ,MBOTQ
      COMMON/PHASES/phi01,phi02,phi0,phiM1,phiM2,phiM3,
     .              phiAT,phiAB,phiATAU,phiAC,phiAS,phiAMU
      COMMON/SFERM3SPEC/MST2,UT,MSB2,UB,MSL2,UTAU,MSNT2
      COMMON/SFERM1SPEC/MSU2,MSD2,MSE2,MSNE2,MSMU2,UMU
      COMMON/CHASPEC/MCH2,U,V
      COMMON/NEUSPEC/MNEU,NEU
      COMMON/HISPEC/MHC,XC,MH0,XH,MA2
      COMMON/CHSFfCOUP/COCHSTbL,COCHSTbR,COCHSBtL,COCHSBtR,COCHSLnL,
     . COCHSNlL,COCHSNlR,COCHSUdL,COCHSUdR,COCHSDuL,COCHSDuR,COCHSEnL,
     . COCHSNeL,COCHSNeR
      COMMON/NEUSFfCOUP/CONESTtL,CONESTtR,CONESBbL,CONESBbR,CONESLlL,
     . CONESLlR,CONESNnL,CONESUuL,CONESUuR,CONESDdL,CONESDdR,CONESEeL,
     . CONESEeR
      COMMON/HINOCOUP/COH0CH,COH0NEU,COHPNEUCHM,COHMNEUCHP


      gg1=dsqrt(g1q)
      gg2=dsqrt(g2q)
      Ytau=mtau/vdq
      DELT(1,1)=1d0
      DELT(1,2)=0d0
      DELT(2,1)=0d0
      DELT(2,2)=1d0


c          A: Inos - Sfermions - fermions

c       I - Charginos - Sfermions - fermions

        DO I=1,2
        DO J=1,2

      COCHSTbL(I,J,1)=Ytq*(V(I,2,1)*UT(J,2,1)+V(I,2,2)*UT(J,2,2))
     .      -gg2*(V(I,1,1)*UT(J,1,1)+V(I,1,2)*UT(J,1,2))
      COCHSTbL(I,J,2)=Ytq*(V(I,2,1)*UT(J,2,2)-V(I,2,2)*UT(J,2,1))
     .      -gg2*(V(I,1,1)*UT(J,1,2)-V(I,1,2)*UT(J,1,1))
      COCHSTbR(I,J,1)=Ybq*(U(I,2,1)*UT(J,1,1)-U(I,2,2)*UT(J,1,2))
      COCHSTbR(I,J,2)=Ybq*(U(I,2,2)*UT(J,1,1)+U(I,2,1)*UT(J,1,2))

      COCHSBtL(I,J,1)=Ybq*(U(I,2,1)*UB(J,2,1)+U(I,2,2)*UB(J,2,2))
     .      -gg2*(U(I,1,1)*UB(J,1,1)+U(I,1,2)*UB(J,1,2))
      COCHSBtL(I,J,2)=Ybq*(U(I,2,1)*UB(J,2,2)-U(I,2,2)*UB(J,2,1))
     .      -gg2*(U(I,1,1)*UB(J,1,2)-U(I,1,2)*UB(J,1,1))
      COCHSBtR(I,J,1)=Ytq*(V(I,2,1)*UB(J,1,1)-V(I,2,2)*UB(J,1,2))
      COCHSBtR(I,J,2)=Ytq*(V(I,2,2)*UB(J,1,1)+V(I,2,1)*UB(J,1,2))

      COCHSLnL(I,J,1)=Ytau*(U(I,2,1)*UTAU(J,2,1)+U(I,2,2)*
     . UTAU(J,2,2))-gg2*(U(I,1,1)*UTAU(J,1,1)+U(I,1,2)*UTAU(J,1,2))
      COCHSLnL(I,J,2)=Ytau*(U(I,2,1)*UTAU(J,2,2)-U(I,2,2)*
     . UTAU(J,2,1))-gg2*(U(I,1,1)*UTAU(J,1,2)-U(I,1,2)*UTAU(J,1,1))

      COCHSUdL(I,J,1)=-gg2*V(I,1,1)*DELT(J,1)
      COCHSUdL(I,J,2)=gg2*V(I,1,2)*DELT(J,1)
      COCHSUdR(I,J,1)=0d0
      COCHSUdR(I,J,2)=0d0

      COCHSDuL(I,J,1)=-gg2*U(I,1,1)*DELT(J,1)
      COCHSDuL(I,J,2)=gg2*U(I,1,2)*DELT(J,1)
      COCHSDuR(I,J,1)=0d0
      COCHSDuR(I,J,2)=0d0

      COCHSEnL(I,J,1)=-gg2*U(I,1,1)*DELT(J,1)
      COCHSEnL(I,J,2)=gg2*U(I,1,2)*DELT(J,1)

        ENDDO

      COCHSNlL(I,1)=-gg2*V(I,1,1)
      COCHSNlL(I,2)=gg2*V(I,1,2)
      COCHSNlR(I,1)=Ytau*U(I,2,1)
      COCHSNlR(I,2)=Ytau*U(I,2,2)

      COCHSNeL(I,1)=-gg2*V(I,1,1)
      COCHSNeL(I,2)=gg2*V(I,1,2)
      COCHSNeR(I,1)=0d0
      COCHSNeR(I,2)=0d0

        ENDDO


c       II - Neutralinos - Sfermions - fermions

        DO I=1,5
        DO J=1,2

      CONESTtL(I,J,1)=-Ytq*(NEU(I,3,1)*UT(J,2,1)+NEU(I,3,2)*UT(J,2,2))
     . -(gg1/3.d0*NEU(I,1,1)+gg2*NEU(I,2,1))/dsqrt(2d0)*UT(J,1,1)
     . -(gg1/3.d0*NEU(I,1,2)+gg2*NEU(I,2,2))/dsqrt(2d0)*UT(J,1,2)
      CONESTtL(I,J,2)=-Ytq*(NEU(I,3,1)*UT(J,2,2)-NEU(I,3,2)*UT(J,2,1))
     . -(gg1/3.d0*NEU(I,1,1)+gg2*NEU(I,2,1))/dsqrt(2d0)*UT(J,1,2)
     . +(gg1/3.d0*NEU(I,1,2)+gg2*NEU(I,2,2))/dsqrt(2d0)*UT(J,1,1)
      CONESTtR(I,J,1)=-Ytq*(NEU(I,3,1)*UT(J,1,1)-NEU(I,3,2)*UT(J,1,2))
     . +2d0*dsqrt(2d0)/3d0*gg1
     .                   *(NEU(I,1,1)*UT(J,2,1)-NEU(I,1,2)*UT(J,2,2))
      CONESTtR(I,J,2)=-Ytq*(NEU(I,3,2)*UT(J,1,1)+NEU(I,3,1)*UT(J,1,2))
     . +2d0*dsqrt(2d0)/3d0*gg1
     .                   *(NEU(I,1,1)*UT(J,2,2)+NEU(I,1,2)*UT(J,2,1))

      CONESBbL(I,J,1)=-Ybq*(NEU(I,4,1)*UB(J,2,1)+NEU(I,4,2)*UB(J,2,2))
     . -(gg1/3.d0*NEU(I,1,1)-gg2*NEU(I,2,1))/dsqrt(2d0)*UB(J,1,1)
     . -(gg1/3.d0*NEU(I,1,2)-gg2*NEU(I,2,2))/dsqrt(2d0)*UB(J,1,2)
      CONESBbL(I,J,2)=-Ybq*(NEU(I,4,1)*UB(J,2,2)-NEU(I,4,2)*UB(J,2,1))
     . -(gg1/3.d0*NEU(I,1,1)-gg2*NEU(I,2,1))/dsqrt(2d0)*UB(J,1,2)
     . +(gg1/3.d0*NEU(I,1,2)-gg2*NEU(I,2,2))/dsqrt(2d0)*UB(J,1,1)
      CONESBbR(I,J,1)=-Ybq*(NEU(I,4,1)*UB(J,1,1)-NEU(I,4,2)*UB(J,1,2))
     .-dsqrt(2d0)/3.d0*gg1*(NEU(I,1,1)*UB(J,2,1)-NEU(I,1,2)*UB(J,2,2))
      CONESBbR(I,J,2)=-Ybq*(NEU(I,4,2)*UB(J,1,1)+NEU(I,4,1)*UB(J,1,2))
     .-dsqrt(2d0)/3.d0*gg1*(NEU(I,1,1)*UB(J,2,2)+NEU(I,1,2)*UB(J,2,1))

      CONESLlL(I,J,1)=-Ytau*(NEU(I,4,1)*UTAU(J,2,1)+NEU(I,4,2)
     . *UTAU(J,2,2))+(gg1*NEU(I,1,1)+gg2*NEU(I,2,1))/dsqrt(2d0)
     . *UTAU(J,1,1)
     . +(gg1*NEU(I,1,2)+gg2*NEU(I,2,2))/dsqrt(2d0)*UTAU(J,1,2)
      CONESLlL(I,J,2)=-Ytau*(NEU(I,4,1)*UTAU(J,2,2)-NEU(I,4,2)
     . *UTAU(J,2,1))+(gg1*NEU(I,1,1)+gg2*NEU(I,2,1))/dsqrt(2d0)
     . *UTAU(J,1,2)
     . -(gg1*NEU(I,1,2)+gg2*NEU(I,2,2))/dsqrt(2d0)*UTAU(J,1,1)
      CONESLlR(I,J,1)=-Ytau*(NEU(I,4,1)*UTAU(J,1,1)-NEU(I,4,2)
     . *UTAU(J,1,2))-dsqrt(2d0)*gg1
     . *(NEU(I,1,1)*UTAU(J,2,1)-NEU(I,1,2)*UTAU(J,2,2))
      CONESLlR(I,J,2)=-Ytau*(NEU(I,4,2)*UTAU(J,1,1)+NEU(I,4,1)
     . *UTAU(J,1,2))-dsqrt(2d0)*gg1
     . *(NEU(I,1,1)*UTAU(J,2,2)+NEU(I,1,2)*UTAU(J,2,1))

      CONESUuL(I,J,1)=-(gg1/3.d0*NEU(I,1,1)+gg2*NEU(I,2,1))
     . /dsqrt(2d0)*DELT(J,1)
      CONESUuL(I,J,2)=(gg1/3.d0*NEU(I,1,2)+gg2*NEU(I,2,2))
     . /dsqrt(2d0)*DELT(J,1)
      CONESUuR(I,J,1)=2d0*dsqrt(2d0)/3d0*gg1*NEU(I,1,1)*DELT(J,2)
      CONESUuR(I,J,2)=2d0*dsqrt(2d0)/3d0*gg1*NEU(I,1,2)*DELT(J,2)

      CONESDdL(I,J,1)=-(gg1/3.d0*NEU(I,1,1)-gg2*NEU(I,2,1))
     . /dsqrt(2d0)*DELT(J,1)
      CONESDdL(I,J,2)=(gg1/3.d0*NEU(I,1,2)-gg2*NEU(I,2,2))
     . /dsqrt(2d0)*DELT(J,1)
      CONESDdR(I,J,1)=-dsqrt(2d0)/3.d0*gg1*NEU(I,1,1)*DELT(J,2)
      CONESDdR(I,J,2)=-dsqrt(2d0)/3.d0*gg1*NEU(I,1,2)*DELT(J,2)

      CONESEeL(I,J,1)=(gg1*NEU(I,1,1)+gg2*NEU(I,2,1))
     . /dsqrt(2d0)*DELT(J,1)
      CONESEeL(I,J,2)=-(gg1*NEU(I,1,2)+gg2*NEU(I,2,2))
     . /dsqrt(2d0)*DELT(J,1)
      CONESEeR(I,J,1)=-dsqrt(2d0)*gg1*NEU(I,1,1)*DELT(J,2)
      CONESEeR(I,J,2)=-dsqrt(2d0)*gg1*NEU(I,1,2)*DELT(J,2)

        ENDDO

      CONESNnL(I,1)=(gg1*NEU(I,1,1)-gg2*NEU(I,2,1))/dsqrt(2d0)
      CONESNnL(I,2)=-(gg1*NEU(I,1,2)-gg2*NEU(I,2,2))/dsqrt(2d0)

        ENDDO


c         B: Inos - Higgs couplings

      DO I=1,5
       XHG(I,1)=XH(I,1)/dsqrt(ZHU)
       XHG(I,2)=XH(I,2)/dsqrt(ZHD)
       XHG(I,3)=XH(I,3)/dsqrt(ZS)
       XHG(I,4)=XH(I,4)*cosb/dsqrt(ZHU)
       XHG(I,5)=XH(I,4)*sinb/dsqrt(ZHD)
       XHG(I,6)=XH(I,5)/dsqrt(ZS)
      ENDDO


c      I - Charginos - Higgs

        DO I=1,5
        DO J=1,2
        DO M=1,2

      COH0CH(I,J,M,1)=(l*DDCOS(Phi01)*(XHG(I,3)*V(J,2,1)*U(M,2,1)
     .       -XHG(I,3)*V(J,2,2)*U(M,2,2)+XHG(I,6)*V(J,2,2)*U(M,2,1)
     .       +XHG(I,6)*V(J,2,1)*U(M,2,2))
     .  +l*DDSIN(Phi01)*(XHG(I,6)*V(J,2,2)*U(M,2,2)
     .       -XHG(I,6)*V(J,2,1)*U(M,2,1)+XHG(I,3)*V(J,2,2)*U(M,2,1)
     .       +XHG(I,3)*V(J,2,1)*U(M,2,2))
     .  +gg2*(XHG(I,1)*(V(J,2,1)*U(M,1,1)-V(J,2,2)*U(M,1,2))
     .       -XHG(I,4)*(V(J,2,2)*U(M,1,1)+V(J,2,1)*U(M,1,2)))
     .  +gg2*(XHG(I,2)*(V(J,1,1)*U(M,2,1)-V(J,1,2)*U(M,2,2))
     .       -XHG(I,5)*(V(J,1,2)*U(M,2,1)+V(J,1,1)*U(M,2,2)))
     .                                                 )/dsqrt(2d0)
      COH0CH(I,J,M,2)=(l*DDCOS(Phi01)*(XHG(I,6)*V(J,2,1)*U(M,2,1)
     .       -XHG(I,6)*V(J,2,2)*U(M,2,2)-XHG(I,3)*V(J,2,2)*U(M,2,1)
     .       -XHG(I,3)*V(J,2,1)*U(M,2,2))
     .  +l*DDSIN(Phi01)*(XHG(I,3)*V(J,2,1)*U(M,2,1)
     .       -XHG(I,3)*V(J,2,2)*U(M,2,2)+XHG(I,5)*V(J,2,2)*U(M,2,1)
     .       +XHG(I,5)*V(J,2,1)*U(M,2,2))
     .  +gg2*(-XHG(I,1)*(V(J,2,1)*U(M,1,2)+V(J,2,2)*U(M,1,1))
     .       -XHG(I,4)*(V(J,2,1)*U(M,1,1)-V(J,2,2)*U(M,1,2)))
     .  +gg2*(-XHG(I,2)*(V(J,1,2)*U(M,2,1)+V(J,1,1)*U(M,2,2))
     .       -XHG(I,5)*(V(J,1,1)*U(M,2,1)-V(J,1,2)*U(M,2,2)))
     .                                                 )/dsqrt(2d0)

        ENDDO
        ENDDO
        ENDDO


c      II - Neutralinos - Higgs

        DO I=1,5
        DO J=1,5
        DO M=1,5

       COH0NEU(I,J,M,1)=(k*DDCOS(Phi02)*
     .    (XHG(I,3)*(NEU(J,5,1)*NEU(M,5,1)-NEU(J,5,2)*NEU(M,5,2))
     .    +XHG(I,6)*(NEU(J,5,1)*NEU(M,5,2)+NEU(J,5,2)*NEU(M,5,1)))
     .  +k*DDSIN(Phi02)*
     .    (-XHG(I,6)*(NEU(J,5,1)*NEU(M,5,1)-NEU(J,5,2)*NEU(M,5,2))
     .    +XHG(I,3)*(NEU(J,5,1)*NEU(M,5,2)+NEU(J,5,2)*NEU(M,5,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,3)*(NEU(J,3,1)*NEU(M,4,1)-NEU(J,3,2)*NEU(M,4,2))
     .    +XHG(I,6)*(NEU(J,3,1)*NEU(M,4,2)+NEU(J,3,2)*NEU(M,4,1)))
     .  -l*DDSIN(Phi01)*
     .    (-XHG(I,6)*(NEU(J,3,1)*NEU(M,4,1)-NEU(J,3,2)*NEU(M,4,2))
     .    +XHG(I,3)*(NEU(J,3,1)*NEU(M,4,2)+NEU(J,3,2)*NEU(M,4,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,1)*(NEU(J,5,1)*NEU(M,4,1)-NEU(J,5,2)*NEU(M,4,2))
     .    +XHG(I,4)*(NEU(J,5,1)*NEU(M,4,2)+NEU(J,5,2)*NEU(M,4,1)))
     .  -l*DDSIN(Phi01)*
     .    (-XHG(I,4)*(NEU(J,5,1)*NEU(M,4,1)-NEU(J,5,2)*NEU(M,4,2))
     .    +XHG(I,1)*(NEU(J,5,1)*NEU(M,4,2)+NEU(J,5,2)*NEU(M,4,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,2)*(NEU(J,5,1)*NEU(M,3,1)-NEU(J,5,2)*NEU(M,3,2))
     .    +XHG(I,5)*(NEU(J,5,1)*NEU(M,3,2)+NEU(J,5,2)*NEU(M,3,1)))
     .  -l*DDSIN(Phi01)*
     .    (-XHG(I,5)*(NEU(J,5,1)*NEU(M,3,1)-NEU(J,5,2)*NEU(M,3,2))
     .    +XHG(I,2)*(NEU(J,5,1)*NEU(M,3,2)+NEU(J,5,2)*NEU(M,3,1)))
     .  +gg1/dsqrt(2d0)*
     .    (XHG(I,1)*(NEU(J,1,1)*NEU(M,3,1)-NEU(J,1,2)*NEU(M,3,2))
     .    -XHG(I,4)*(NEU(J,1,1)*NEU(M,3,2)+NEU(J,1,2)*NEU(M,3,1))
     .    -XHG(I,2)*(NEU(J,1,1)*NEU(M,4,1)-NEU(J,1,2)*NEU(M,4,2))
     .    +XHG(I,5)*(NEU(J,1,1)*NEU(M,4,2)+NEU(J,1,2)*NEU(M,4,1)))
     .  -gg2/dsqrt(2d0)*
     .    (XHG(I,1)*(NEU(J,2,1)*NEU(M,3,1)-NEU(J,2,2)*NEU(M,3,2))
     .    -XHG(I,4)*(NEU(J,2,1)*NEU(M,3,2)+NEU(J,2,2)*NEU(M,3,1))
     .    -XHG(I,2)*(NEU(J,2,1)*NEU(M,4,1)-NEU(J,2,2)*NEU(M,4,2))
     .    +XHG(I,5)*(NEU(J,2,1)*NEU(M,4,2)+NEU(J,2,2)*NEU(M,4,1)))
     .                                               )/dsqrt(2d0)

       COH0NEU(I,J,M,2)=(k*DDSIN(Phi02)*
     .    (XHG(I,3)*(NEU(J,5,1)*NEU(M,5,1)-NEU(J,5,2)*NEU(M,5,2))
     .    +XHG(I,6)*(NEU(J,5,1)*NEU(M,5,2)+NEU(J,5,2)*NEU(M,5,1)))
     .  +k*DDCOS(Phi02)*
     .    (XHG(I,6)*(NEU(J,5,1)*NEU(M,5,1)-NEU(J,5,2)*NEU(M,5,2))
     .    -XHG(I,3)*(NEU(J,5,1)*NEU(M,5,2)+NEU(J,5,2)*NEU(M,5,1)))
     .  -l*DDSIN(Phi01)*
     .    (XHG(I,3)*(NEU(J,3,1)*NEU(M,4,1)-NEU(J,3,2)*NEU(M,4,2))
     .    +XHG(I,6)*(NEU(J,3,1)*NEU(M,4,2)+NEU(J,3,2)*NEU(M,4,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,6)*(NEU(J,3,1)*NEU(M,4,1)-NEU(J,3,2)*NEU(M,4,2))
     .    -XHG(I,3)*(NEU(J,3,1)*NEU(M,4,2)+NEU(J,3,2)*NEU(M,4,1)))
     .  -l*DDSIN(Phi01)*
     .    (XHG(I,1)*(NEU(J,5,1)*NEU(M,4,1)-NEU(J,5,2)*NEU(M,4,2))
     .    +XHG(I,4)*(NEU(J,5,1)*NEU(M,4,2)+NEU(J,5,2)*NEU(M,4,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,4)*(NEU(J,5,1)*NEU(M,4,1)-NEU(J,5,2)*NEU(M,4,2))
     .    -XHG(I,1)*(NEU(J,5,1)*NEU(M,4,2)+NEU(J,5,2)*NEU(M,4,1)))
     .  -l*DDSIN(Phi01)*
     .    (XHG(I,2)*(NEU(J,5,1)*NEU(M,3,1)-NEU(J,5,2)*NEU(M,3,2))
     .    +XHG(I,5)*(NEU(J,5,1)*NEU(M,3,2)+NEU(J,5,2)*NEU(M,3,1)))
     .  -l*DDCOS(Phi01)*
     .    (XHG(I,5)*(NEU(J,5,1)*NEU(M,3,1)-NEU(J,5,2)*NEU(M,3,2))
     .    -XHG(I,2)*(NEU(J,5,1)*NEU(M,3,2)+NEU(J,5,2)*NEU(M,3,1)))
     .  -gg1/dsqrt(2d0)*
     .    (XHG(I,1)*(NEU(J,1,1)*NEU(M,3,2)+NEU(J,1,2)*NEU(M,3,1))
     .    +XHG(I,4)*(NEU(J,1,1)*NEU(M,3,1)-NEU(J,1,2)*NEU(M,3,2))
     .    -XHG(I,2)*(NEU(J,1,1)*NEU(M,4,2)+NEU(J,1,2)*NEU(M,4,1))
     .    -XHG(I,5)*(NEU(J,1,1)*NEU(M,4,1)-NEU(J,1,2)*NEU(M,4,2)))
     .  +gg2/dsqrt(2d0)*
     .    (XHG(I,1)*(NEU(J,2,1)*NEU(M,3,2)+NEU(J,2,2)*NEU(M,3,1))
     .    +XHG(I,4)*(NEU(J,2,1)*NEU(M,3,1)-NEU(J,2,2)*NEU(M,3,2))
     .    -XHG(I,2)*(NEU(J,2,1)*NEU(M,4,2)+NEU(J,2,2)*NEU(M,4,1))
     .    -XHG(I,5)*(NEU(J,2,1)*NEU(M,4,1)-NEU(J,2,2)*NEU(M,4,2)))
     .                                               )/dsqrt(2d0)

        ENDDO
        ENDDO

        DO J=1,5
        DO M=J,5

      COH0NEU(I,J,M,1)=COH0NEU(I,J,M,1)+COH0NEU(I,M,J,1)
      COH0NEU(I,M,J,1)=COH0NEU(I,J,M,1)
      COH0NEU(I,J,M,2)=COH0NEU(I,J,M,2)+COH0NEU(I,M,J,2)
      COH0NEU(I,M,J,2)=COH0NEU(I,J,M,2)

        ENDDO
        ENDDO
        ENDDO


c      III - Neutralinos - Charginos - charged Higgs

        DO I=1,2
        DO J=1,5
        DO M=1,2

      COHPNEUCHM(I,J,M,1)=l*XC(I,1)*
     .     (DDCOS(Phi01)*(NEU(J,5,1)*U(M,2,1)-NEU(J,5,2)*U(M,2,2))
     .     +DDSIN(Phi01)*(NEU(J,5,1)*U(M,2,2)+NEU(J,5,2)*U(M,2,1)))
     . -(gg1*NEU(J,1,1)+gg2*NEU(J,2,1))/dsqrt(2d0)*XC(I,2)*U(M,2,1)
     . +(gg1*NEU(J,1,2)+gg2*NEU(J,2,2))/dsqrt(2d0)*XC(I,2)*U(M,2,2)
     . +gg2*XC(I,2)*(NEU(J,4,1)*U(M,1,1)-NEU(J,4,2)*U(M,1,2))
      COHPNEUCHM(I,J,M,2)=l*XC(I,1)*
     .     (DDSIN(Phi01)*(NEU(J,5,1)*U(M,2,1)-NEU(J,5,2)*U(M,2,2))
     .     -DDCOS(Phi01)*(NEU(J,5,1)*U(M,2,2)+NEU(J,5,2)*U(M,2,1)))
     . +(gg1*NEU(J,1,1)+gg2*NEU(J,2,1))*XC(I,2)/dsqrt(2d0)*U(M,2,2)
     . +(gg1*NEU(J,1,2)+gg2*NEU(J,2,2))*XC(I,2)/dsqrt(2d0)*U(M,2,1)
     . -gg2*XC(I,2)*(NEU(J,4,1)*U(M,1,2)+NEU(J,4,2)*U(M,1,1))
      COHMNEUCHP(I,J,M,1)=l*XC(I,2)*
     .      (DDCOS(Phi01)*(NEU(J,5,1)*V(M,2,1)-NEU(J,5,2)*V(M,2,2))
     .      +DDSIN(Phi01)*(NEU(J,5,1)*V(M,2,2)+NEU(J,5,2)*V(M,2,1)))
     . +(gg1*NEU(J,1,1)+gg2*NEU(J,2,1))*XC(I,1)/dsqrt(2d0)*V(M,2,1)
     . -(gg1*NEU(J,1,2)+gg2*NEU(J,2,2))*XC(I,1)/dsqrt(2d0)*V(M,2,2)
     . +gg2*XC(I,1)*(NEU(J,3,1)*V(M,1,1)-NEU(J,3,2)*V(M,1,2))
      COHMNEUCHP(I,J,M,2)=l*XC(I,2)*
     .      (DDSIN(Phi01)*(NEU(J,5,1)*V(M,2,1)-NEU(J,5,2)*V(M,2,2))
     .      -DDCOS(Phi01)*(NEU(J,5,1)*V(M,2,2)+NEU(J,5,2)*V(M,2,1)))
     . -(gg1*NEU(J,1,1)+gg2*NEU(J,2,1))*XC(I,1)/dsqrt(2d0)*V(M,2,2)
     . -(gg1*NEU(J,1,2)+gg2*NEU(J,2,2))*XC(I,1)/dsqrt(2d0)*V(M,2,1)
     . -gg2*XC(I,1)*(NEU(J,3,1)*V(M,1,2)+NEU(J,3,2)*V(M,1,1))

        ENDDO
        ENDDO
        ENDDO


c         C: Sfermions - Higgs couplings

c These have already been computed within mhiggsloop_pole_CPV.f


      RETURN
      END
