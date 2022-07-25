*   Subroutines to integrate the RGEs for the gauge and Yukawa couplings
*   above the messenger scale: the 1-loop coefficients of the gauge beta
*   functions include the contributions from n5 pairs of messengers in
*   5- and 5_bar-representations of SU(5)

      SUBROUTINE
     .   ODEINTGMGUT(YSTART,NVAR,X1,X2,EPS,DERIVSGMGUT,RKQSGMGUT,IFAIL)

*   Driver subroutine to integrate Ordinary Differential Equations
*   using the Runge-Kutta 5th order method with adaptative step size
*   (from Numerical Recipes)
*   IFAIL=1 stepsize smaller than minimum in ODEINTGMGUT
*   IFAIL=2 too many steps in ODEINTGMGUT
*   IFAIL=3 stepsize underflow in RKQSGMGUT
*   IFAIL=4 max(|dy(i)/dx|*(1+x)/(1+|y(i)|)) > 1/eps^2

      IMPLICIT NONE

      INTEGER NVAR,NMAX,IFAIL,MAXSTP,NSTP,I
      PARAMETER(MAXSTP=10000,NMAX=500)

      DOUBLE PRECISION YSTART(NVAR),Y(NMAX),YSCAL(NMAX)
      DOUBLE PRECISION DYDX(NMAX),X,X1,X2,EPS,TINY
      DOUBLE PRECISION H1,HMIN,H,HDID,HNEXT

      EXTERNAL DERIVSGMGUT,RKQSGMGUT

      IFAIL=0
      TINY=EPS**4
      H1=DSQRT(EPS)*DABS(X2-X1)
      HMIN=EPS**2*DABS(X2-X1)

      X=X1
      H=SIGN(H1,X2-X1)
      DO I=1,NVAR
       Y(I)=YSTART(I)
      ENDDO

      DO NSTP=1,MAXSTP

       CALL DERIVSGMGUT(NVAR,X,Y,DYDX)

       DO I=1,NVAR
        IF(DABS(DYDX(I))*(1d0+X)/(1d0+DABS(Y(I))).GT.EPS**(-2))IFAIL=4
       ENDDO
       IF(IFAIL.GT.0)RETURN

       DO I=1,NVAR
        YSCAL(I)=DABS(Y(I))+DABS(H*DYDX(I))+TINY
       ENDDO

       CALL RKQSGMGUT(Y,DYDX,NVAR,X,X2,H,EPS,YSCAL,
     .  HDID,HNEXT,DERIVSGMGUT,IFAIL)
       IF(IFAIL.GT.0)RETURN
      
       IF(DABS(Y(2)-5d0/3d0*Y(1)).LT.EPS)THEN
        X2=X
        DO I=1,NVAR
         YSTART(I)=Y(I)
        ENDDO
        RETURN
       ENDIF

       IF(ABS(HNEXT).LT.HMIN)THEN
        X2=X
        DO I=1,NVAR
         YSTART(I)=Y(I)
        ENDDO
        IFAIL=1
        RETURN
       ENDIF

       H=HNEXT

      ENDDO

      IFAIL=2

      RETURN
      END


      SUBROUTINE RKQSGMGUT(Y,DYDX,N,X,X2,HTRY,EPS,YSCAL,HDID,
     .    HNEXT,DERIVSGMGUT,IFAIL)

*   Stepper subroutine for ODEINTGMGUT

      IMPLICIT NONE

      INTEGER IFAIL,I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION EPS,HDID,HNEXT,HTRY,X,X2,DYDX(N),Y(N),YSCAL(N)
      DOUBLE PRECISION ERRMAX,H,HTEMP,XNEW,YERR(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION SAFETY,PGROW,PSHRNK,ERRCON,G,GTEMP

      EXTERNAL DERIVSGMGUT

      SAFETY=.9d0
      PGROW=-.2d0
      PSHRNK=-.25d0
      ERRCON=(5d0/SAFETY)**(1d0/PGROW)
      H=HTRY
      G=Y(2)-5d0/3d0*Y(1)

1     CALL RKCKGMGUT(Y,DYDX,N,X,H,YTEMP,YERR,DERIVSGMGUT)

      ERRMAX=0d0
      DO I=1,N
       ERRMAX=MAX(ERRMAX,DABS(YERR(I)/YSCAL(I)))
      ENDDO
      ERRMAX=ERRMAX/EPS

      IF(ERRMAX.GT.1d0)THEN
       HTEMP=SAFETY*H*(ERRMAX**PSHRNK)
       H=SIGN(MAX(DABS(HTEMP),.1d0*DABS(H)),H)
       XNEW=X+H
       IF(XNEW.EQ.X)THEN
        IFAIL=3
        RETURN
       ENDIF
       GOTO 1
      ENDIF

      GTEMP=YTEMP(2)-5d0/3d0*YTEMP(1)
      IF(GTEMP.LT.-EPS)THEN
       IFAIL=-1
       X2=X+H
       H=G/(G-GTEMP)*H
       GOTO 1
      ENDIF

      HDID=H
      X=X+H
      DO I=1,N
       Y(I)=YTEMP(I)
      ENDDO

      IF(ERRMAX.GT.ERRCON)THEN
       HNEXT=SAFETY*H*(ERRMAX**PGROW)
      ELSE
        HNEXT=5d0*H
      ENDIF
      
      IF(X+HNEXT.GT.X2 .AND. IFAIL.EQ.-1) HNEXT=X2-X

      RETURN

      END


      SUBROUTINE RKCKGMGUT(Y,DYDX,N,X,H,YOUT,YERR,DERIVSGMGUT)

*   Algorithm subroutine for ODEINTGMGUT

      IMPLICIT NONE

      INTEGER I,N,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION H,X,DYDX(N),Y(N),YERR(N),YOUT(N)
      DOUBLE PRECISION AK2(NMAX),AK3(NMAX),AK4(NMAX)
      DOUBLE PRECISION AK5(NMAX),AK6(NMAX),YTEMP(NMAX)
      DOUBLE PRECISION A2,A3,A4,A5,A6,C1,C3,C4,C6
      DOUBLE PRECISION B21,B31,B32,B41,B42,B43,B51,B52,B53,B54
      DOUBLE PRECISION B61,B62,B63,B64,B65,DC1,DC3,DC4,DC5,DC6

      PARAMETER(A2=.2d0, A3=.3d0, A4=.6d0, A5=1d0, A6=.875d0,
     .       B21=.2d0, B31=3d0/40d0, B32=9d0/40d0, B41=.3d0, B42=-.9d0,
     .       B43=1.2d0, B51=-11d0/54d0, B52=2.5d0, B53=-70d0/27d0,
     .       B54=35d0/27d0, B61=1631d0/55296d0, B62=175d0/512d0,
     .       B63=575d0/13824d0, B64=44275d0/110592d0,
     .       B65=253d0/4096d0, C1=37d0/378d0, C3=250d0/621d0,
     .       C4=125d0/594d0, C6=512d0/1771d0, DC1=C1-2825d0/27648d0,
     .       DC3=C3-18575d0/48384d0, DC4=C4-13525d0/55296d0,
     .       DC5=-277d0/14336d0, DC6=C6-.25d0)

      EXTERNAL DERIVSGMGUT

      DO I=1,N
       YTEMP(I)=Y(I)+B21*H*DYDX(I)
      ENDDO
      CALL DERIVSGMGUT(N,X+A2*H,YTEMP,AK2)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B31*DYDX(I)+B32*AK2(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A3*H,YTEMP,AK3)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B41*DYDX(I)+B42*AK2(I)+B43*AK3(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A4*H,YTEMP,AK4)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B51*DYDX(I)+B52*AK2(I)+B53*AK3(I)+B54*AK4(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A5*H,YTEMP,AK5)
      DO I=1,N
       YTEMP(I)=Y(I)+H*(B61*DYDX(I)+B62*AK2(I)+B63*AK3(I)+B64*AK4(I)+
     .  B65*AK5(I))
      ENDDO
      CALL DERIVSGMGUT(N,X+A6*H,YTEMP,AK6)
      DO I=1,N
       YOUT(I)=Y(I)+H*(C1*DYDX(I)+C3*AK3(I)+C4*AK4(I)+C6*AK6(I))
      ENDDO
      DO I=1,N
       YERR(I)=H*(DC1*DYDX(I)+DC3*AK3(I)+DC4*AK4(I)+DC5*AK5(I)+DC6*
     .  AK6(I))
      ENDDO

      RETURN
      END


      SUBROUTINE DERIVSGMGUT(N,X,Y,F)

*   2-loop Renormalization group equations for G1, G2, G3,
*   lambda, kappa, htop, hbot, htau, lambda_PhiPhi, lambda_TT,
*   lambda_Hu, lambda_Hd, lambda_top, lambda_bot, lambda_tau
*   to be integrated by ODEINTGMGUT

      IMPLICIT NONE

      INTEGER N
      DOUBLE PRECISION X,Y(N),F(N),PI,c2
      DOUBLE PRECISION G1,G2,G3,L,K,HT,HB,HL
      DOUBLE PRECISION LPP,LTT,LU,LD,LT,LB,LL
      DOUBLE PRECISION DL,DK,DHT,DHB,DHL,DLPP
      DOUBLE PRECISION DLTT,DLU,DLD,DLT,DLB,DLL
      DOUBLE PRECISION MSUSYEFF,MMESS,N5

      COMMON/MESCAL/MSUSYEFF,MMESS,N5

      PI=4d0*DATAN(1d0)
      c2=1d0/(16d0*PI**2)

      G1=Y(1)
      G2=Y(2)
      G3=Y(3)
      L=Y(4)**2
      K=Y(5)**2
      HT=Y(6)**2
      HB=Y(7)**2
      HL=Y(8)**2
      LPP=Y(9)**2
      LTT=Y(10)**2
      LU=Y(11)**2
      LD=Y(12)**2
      LT=Y(13)**2
      LB=Y(14)**2
      LL=Y(15)**2
      DL=Y(4)
      DK=Y(5)
      DHT=Y(6)
      DHB=Y(7)
      DHL=Y(8)
      DLPP=Y(9)
      DLTT=Y(10)
      DLU=Y(11)
      DLD=Y(12)
      DLT=Y(13)
      DLB=Y(14)
      DLL=Y(15)

      F(1)= (N5*5d0/3d0+11d0)*G1**2
     .    + c2*G1**2*(269d0/9d0*G1 + 27d0*G2 + 152d0/3d0*G3
     .    - 26d0/3d0*HT - 14d0/3d0*HB- 6d0*HL - 2d0*L
     .    - 4d0/3d0*LTT - 2d0*LPP - 2d0*LU - 2d0*LD
     .    - 26d0/3d0*LT - 14d0/3d0*LB - 6*LL
     .    - N5*(140d0/27d0*G1 + 12d0*G2 + 128d0/9d0*G3)
     .    + N5**2*(35d0/27d0*G1 + 3d0*G2 + 32d0/9d0*G3))

      F(2)= (N5+1d0)*G2**2
     .    + c2*G2**2*( 9d0*G1 + 43d0*G2 + 24d0*G3
     .    - 6d0*HT - 6d0*HB - 2d0*HL - 2d0*L
     .    - 2d0*LPP - 2d0*LU - 2d0*LD - 6d0*LT - 6d0*LB - 2d0*LL
     .    - N5*(4d0*G1 + 8d0*G2) + N5**2*(G1 + 3d0*G2))

      F(3)= (N5-3d0)*G3**2
     .    + c2*G3**2*(19d0/3d0*G1 + 9d0*G2 + 46d0*G3
     .    - 4d0*HT - 4d0*HB - 2d0*LTT - 4d0*LT - 4d0*LB
     .    - N5*(16d0/9d0*G1 + 46d0/3d0*G3)
     .    + N5**2*(4d0/9d0*G1 + 16d0/3d0*G3))

      F(4)= (DL*(-G1 - 3d0*G2 + 4d0*L + 2d0*K + 3d0*HT
     .    + 3d0*HB + HL + 2d0*LPP + 3d0*LTT + 4d0*LU + 4d0*LD)
     .    + 3d0*DLD*DLT*DHT + 3d0*DLU*DLB*DHB
     .    + DLU*DLL*DHL
     .    + c2*(DL*(23d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .    - 2d0/3d0*G1*HB + 16d0*G3*HB - 9d0*HB**2 + 2d0*G1*HL
     .    - 3d0*HL**2 + 4d0/3d0*G1*HT + 16d0*G3*HT - 6d0*HB*HT
     .    - 9d0*HT**2- 8d0*K**2 + 2d0*G1*L + 6d0*G2*L - 9d0*HB*L
     .    - 3d0*HL*L - 9d0*HT*L - 12d0*K*L - 10d0*L**2
     .    - 9d0*HB*LB - 3d0*HT*LB + 2d0*G1*LD + 6d0*G2*LD
     .    - 9d0*HB*LD - 3d0*HL*LD - 12d0*K*LD - 20d0*L*LD
     .    - 10d0*LD**2 - 3d0*HL*LL + 2d0*G1*LPP + 6d0*G2*LPP
     .    - 8d0*K*LPP - 4d0*L*LPP - 4d0*LD*LPP - 4d0*LPP**2
     .    - 3d0*HB*LT - 9d0*HT*LT - 9d0*LD*LT + 4d0/3d0*G1*LTT
     .    + 16d0*G3*LTT - 12d0*K*LTT - 6d0*L*LTT - 6d0*LD*LTT
     .    - 6d0*LTT**2 + 2d0*G1*LU + 6d0*G2*LU - 9d0*HT*LU
     .    - 12d0*K*LU - 20d0*L*LU - 9d0*LB*LU - 10d0*LD*LU
     .    - 3d0*LL*LU - 4d0*LPP*LU - 6d0*LTT*LU - 10d0*LU**2
     .    + N5*(5d0/3d0*G1**2 + 3d0*G2**2))
     .    + DHT*DLT*DLD*(4d0/3d0*G1 + 16d0*G3 - 9d0*HT - 3d0*HB
     .    - 18d0*L - 3d0*LU - 9d0*LT - 3d0*LB)
     .    + DHB*DLB*DLU*(-2d0/3d0*G1 + 16d0*G3 - 3d0*HT - 9d0*HB
     .    - 18d0*L  - 3d0*LD - 3d0*LT - 9d0*LB)
     .    + DHL*DLL*DLU*(2d0*G1 - 3d0*HL - 6d0*L - LD
     .    - 3d0*LL)))/2d0

      F(5)= 3d0*DK*(L + K + LPP + 1.5d0*LTT + LU + LD)
     .    + c2*DK*(-24d0*K**2 + 6d0*G1*L + 18d0*G2*L - 18d0*HB*L
     .    - 6d0*HL*L - 18*HT*L - 24*K*L + 6*G1*LD - 24d0*L*LD
     .    + 18d0*G2*LD- 18d0*HB*LD - 6d0*HL*LD - 24d0*K*LD - 12*L**2
     .    - 12d0*LD**2 + 6d0*G1*LPP + 18d0*G2*LPP - 24d0*K*LPP
     .    - 12d0*LPP**2 - 18d0*LD*LT + 4d0*G1*LTT + 48d0*G3*LTT
     .    - 36d0*K*LTT - 18d0*LTT**2 + 6d0*G1*LU + 18d0*G2*LU
     .    - 18d0*HT*LU - 24d0*K*LU - 24d0*L*LU - 18d0*LB*LU
     .    - 6d0*LL*LU - 12d0*LU**2 - 36d0*DHT*DL*DLD*DLT
     .    - 36d0*DHB*DL*DLB*DLU - 12d0*DHL*DL*DLL*DLU)/2d0

      F(6)= (DHT*(-13d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .    + L + 6d0*HT + HB + LU + 6d0*LT + LB)
     .    + DL*DLD*DLT
     .    + c2*(DHT*(2743d0/162d0*G1**2 + 5d0/3d0*G1*G2 + 15d0/2d0*G2**2
     .    + 136d0/27d0*G1*G3 + 8d0*G2*G3 - 16d0/9d0*G3**2 + 2d0*G1*HT
     .     + 2d0/3d0*G1*HB- 5d0*HB**2 - HB*HL + 6d0*G2*HT + 16d0*G3*HT
     .    - 5d0*HB*HT - 22d0*HT**2 - 4d0*HB*L - HL*L - 3d0*HT*L
     .    - 2d0*K*L - 3d0*L**2 + 2d0/3d0*G1*LB - 10d0*HB*LB - 5d0*HT*LB
     .    - 5d0*LB**2 - HB*LD - 3d0*L*LD - LB*LL - 2d0*L*LPP
     .    + 2d0*G1*LT + 6d0*G2*LT + 16d0*G3*LT - 5d0*HB*LT - 44d0*HT*LT
     .    - 5d0*LB*LT - 3d0*LD*LT - 22d0*LT**2 - 3d0*L*LTT
     .    - 2d0*DHB*DHL*DLB*DLL - 6d0*DHT*DL*DLD*DLT
     .    - 8d0*DHB*DL*DLB*DLU - 2d0*DHL*DL*DLL*DLU
     .    - 3d0*HT*LU - 2d0*K*LU - 6d0*L*LU - 4d0*LB*LU - 2d0*LD*LU
     .    - LL*LU - 2d0*LPP*LU - 3d0*LTT*LU - 3d0*LU**2
     .    + N5*(65d0/27d0*G1**2 + 3d0*G2**2 + 16d0/3d0*G3**2))
     .    - DL*DLD*DLT*(3d0*HB + HL + 3d0*L + 2d0*K
     .    + 2d0*LPP + 3d0*LTT + 3d0*LU + 3d0*LD)
     .    - 3d0*DHB*DLB*DLD*DLT*DLU - DHL*DLD*DLL*DLT*DLU))/2d0

      F(7)= (DHB*(-7d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .    + L + HT + 6d0*HB + HL + LD + 6d0*LB + LT)
     .    + DL*DLU*DLB + DLB*DLL*DHL
     .    + c2*(DHB*(1435d0/162d0*G1**2 + 5d0/3d0*G1*G2 + 15d0/2d0*G2**2
     .    + 40d0/27d0*G1*G3 + 8d0*G2*G3 - 16d0/9d0*G3**2 + 2d0/3d0*G1*HB
     .    + 6d0*G2*HB + 16d0*G3*HB - 22d0*HB**2 + 2d0*G1*HL
     .    - 3d0*HB*HL - 3d0*HL**2 + 4d0/3d0*G1*HT - 5d0*HB*HT
     .    - 5d0*HT**2 - 3d0*HB*L - 4d0*HT*L - 2d0*K*L - 3d0*L**2
     .    + 2d0/3d0*G1*LB + 6d0*G2*LB + 16d0*G3*LB - 44d0*HB*LB
     .    - 5d0*HT*LB - 22d0*LB**2 - 3d0*HB*LD - 2d0*K*LD
     .    - 6d0*L*LD - 3d0*LD**2 - 3d0*HL*LL - 3d0*LB*LL
     .    - 2d0*L*LPP - 2d0*LD*LPP + 4d0/3d0*G1*LT - 5d0*HB*LT
     .    - 10d0*HT*LT - 5d0*LB*LT - 4d0*LD*LT - 5d0*LT**2 - 3d0*L*LTT
     .    - 3d0*LD*LTT - HT*LU - 3d0*L*LU - 3d0*LB*LU - 2d0*LD*LU
     .    + N5*(35d0/27d0*G1**2 + 3d0*G2**2 + 16d0/3d0*G3**2)
     .    - 6d0*DHB*DHL*DLB*DLL - 8d0*DL*DHT*DLD*DLT
     .    - 6d0*DL*DHB*DLU*DLB)
     .    - DHL*DLB*DLL*(-2d0*G1 + 3d0*HL + 3d0*LL)
     .    - DL*DLB*DLU*(3d0*HT + 3d0*L + 2d0*K + 2d0*LPP + 3d0*LTT
     .    + 3d0*LU + 3d0*LD) - 3d0*DHT*DLB*DLD*DLT*DLU))/2d0

      F(8)= (DHL*(-3d0*G1 - 3d0*G2
     .    + L + 3d0*HB + 4d0*HL + LD + 4*LL)
     .    + DL*DLU*DLL + 3d0*DLB*DLL*DHB
     .    + c2*(DHL*(75d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .    - 2d0/3d0*G1*HB + 16d0*G3*HB - 9d0*HB**2 + 2d0*G1*HL
     .    + 6d0*G2*HL - 9d0*HB*HL - 10d0*HL**2 - 3d0*HB*HT - 3d0*HL*L
     .    - 3d0*HT*L - 2d0*K*L - 3d0*L**2 - 9d0*HB*LB - 3d0*HL*LD
     .    - 2d0*K*LD - 6d0*L*LD - 3d0*LD**2 + 2d0*G1*LL + 6d0*G2*LL
     .    - 20d0*HL*LL - 9d0*LB*LL - 10d0*LL**2 - 2d0*L*LPP
     .    - 2d0*LD*LPP - 3d0*HB*LT - 3d0*LD*LT - 3d0*L*LTT
     .    - 3d0*LD*LTT - 3d0*L*LU - 2d0*LD*LU - 3d0*LL*LU
     .    + N5*(5d0*G1**2 + 3d0*G2**2) - 18d0*DHB*DHL*DLB*DLL
     .    - 6d0*DHT*DL*DLD*DLT - 6d0*DHL*DL*DLL*DLU)
     .    - DHB*DLB*DLL*(2d0/3d0*G1 - 16d0*G3 + 3d0*HT + 9d0*HB
     .    + 3d0*LT + 9d0*LB)
     .    - DL*DLL*DLU*(3d0*HT + 2d0*K + 3d0*L + 3d0*LU + 3d0*LD
     .    + 2d0*LPP + 3d0*LTT)- 3d0*DHT*DLD*DLL*DLT*DLU))/2d0

      F(9)= DLPP*(-G1 - 3d0*G2 + 2d0*L + 2d0*K
     .    + 4d0*LPP + 3d0*LTT + 2d0*LU + 2d0*LD)/2d0
     .    + c2*DLPP*(23d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .    - 8d0*K**2 + 2d0*G1*L + 6d0*G2*L - 6d0*HB*L - 2d0*HL*L
     .    - 6d0*HT*L - 8d0*K*L - 4d0*L**2 + 2d0*G1*LD + 6d0*G2*LD
     .    - 6d0*HB*LD - 2d0*HL*LD - 8d0*K*LD - 8d0*L*LD - 4d0*LD**2
     .    + 2d0*G1*LPP + 6d0*G2*LPP - 12d0*K*LPP - 4d0*L*LPP
     .    - 4d0*LD*LPP - 10d0*LPP**2 - 6d0*LD*LT + 4d0/3d0*G1*LTT
     .    + 16d0*G3*LTT - 12d0*K*LTT - 6d0*LPP*LTT - 6d0*LTT**2
     .    + 2d0*G1*LU + 6d0*G2*LU - 6d0*HT*LU - 8d0*K*LU - 8d0*L*LU
     .    - 6d0*LB*LU - 2d0*LL*LU - 4d0*LPP*LU - 4d0*LU**2
     .    + N5*(5d0/3d0*G1**2 + 3d0*G2**2) - 12d0*DHT*DL*DLD*DLT
     .    - 12d0*DHB*DL*DLB*DLU - 4d0*DHL*DL*DLL*DLU)/2d0

      F(10)= DLTT*(-4d0/9d0*G1 - 16d0/3d0*G3 + 2d0*L + 2d0*K
     .     + 2d0*LPP + 5d0*LTT + 2d0*LU + 2d0*LD)/2d0
     .     + c2*DLTT*(404d0/81d0*G1**2 + 64d0/27d0*G1*G3 - 8d0*K**2
     .     - 16d0/9d0*G3**2 + 2d0*G1*L + 6d0*G2*L - 6d0*HB*L - 2d0*HL*L
     .     - 6d0*HT*L - 8d0*K*L - 4d0*L**2 + 2d0*G1*LD + 6d0*G2*LD
     .     - 6d0*HB*LD - 2d0*HL*LD - 8d0*K*LD - 8d0*L*LD - 4d0*LD**2
     .     + 2d0*G1*LPP + 6d0*G2*LPP - 8d0*K*LPP - 4d0*LPP**2
     .     - 6d0*LD*LT + 4d0/3d0*G1*LTT + 16d0*G3*LTT - 16d0*K*LTT
     .     - 4d0*L*LTT - 4d0*LD*LTT - 4d0*LPP*LTT - 14d0*LTT**2
     .     + 2d0*G1*LU + 6d0*G2*LU - 6d0*HT*LU - 8d0*K*LU - 8d0*L*LU
     .     - 6d0*LB*LU - 2d0*LL*LU - 4d0*LTT*LU - 4d0*LU**2
     .     + N5*(20d0/27d0*G1**2 + 16d0/3d0*G3**2)
     .     - 12d0*DHB*DL*DLB*DLU - 4d0*DHL*DL*DLL*DLU
     .     - 12d0*DHT*DL*DLD*DLT)/2d0

      F(11)= (DLU*(-G1 - 3d0*G2 + 4d0*L + 2d0*K
     .    + 2d0*LPP + 3d0*LTT + 4d0*LU + 2d0*LD
     .    + 3d0*HT + 3d0*LB + LL)
     .    + 3d0*DL*DLB*DHB + DL*DLL*DHL
     .    + c2*(DLU*(23d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .    + 4d0/3d0*G1*HT + 16d0*G3*HT - 3d0*HB*HT - 9d0*HT**2
     .    - 8d0*K**2 + 2d0*G1*L + 6d0*G2*L - 9d0*HB*L - 3d0*HL*L
     .    - 9d0*HT*L - 12d0*K*L - 10d0*L**2 - 2d0/3d0*G1*LB
     .    + 16d0*G3*LB - 9d0*HB*LB - 6d0*HT*LB - 9d0*LB**2
     .    + 2d0*G1*LD + 6d0*G2*LD - 6d0*HB*LD - 2d0*HL*LD
     .    - 8d0*K*LD - 14d0*L*LD - 4d0*LD**2 + 2d0*G1*LL
     .    - 3d0*HL*LL - 3d0*LL**2 + 2d0*G1*LPP + 6d0*G2*LPP
     .    - 8d0*K*LPP - 4d0*L*LPP - 4d0*LPP**2 - 9d0*HT*LT
     .    - 3d0*LB*LT - 6d0*LD*LT + 4d0/3d0*G1*LTT + 16d0*G3*LTT
     .    - 12d0*K*LTT - 6d0*L*LTT - 6d0*LTT**2 + 2d0*G1*LU
     .    + 6d0*G2*LU - 9d0*HT*LU - 12d0*K*LU - 20d0*L*LU
     .    - 9d0*LB*LU - 4d0*LD*LU - 3d0*LL*LU - 4d0*LPP*LU
     .    - 6d0*LTT*LU - 10d0*LU**2 + N5*(5d0/3d0*G1**2 + 3d0*G2**2)
     .    - 18d0*DHB*DL*DLB*DLU - 6d0*DHL*DL*DLL*DLU
     .    - 15d0*DHT*DL*DLD*DLT)
     .    - DHB*DL*DLB*(2d0/3d0*G1 - 16d0*G3 + 9d0*HB + 3d0*HT
     .    + 9d0*LB + 3d0*LT)
     .    + DL*DHL*DLL*(2d0*G1 - 3d0*HL - 3d0*LL)))/2d0

      F(12)= (DLD*(-G1 - 3d0*G2 + 4d0*L + 2d0*K
     .    + 2d0*LPP + 3d0*LTT + 2d0*LU + 4d0*LD
     .    + 3d0*HB + HL + 3d0*LT)
     .    + 3d0*DL*DLT*DHT
     .    + c2*(DLD*(23d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .    - 2d0/3d0*G1*HB + 16d0*G3*HB - 9d0*HB**2 + 2d0*G1*HL
     .    - 3d0*HL**2 - 3d0*HB*HT - 8d0*K**2 + 2d0*G1*L + 6d0*G2*L
     .    - 9d0*HB*L - 3d0*HL*L - 9d0*HT*L - 12d0*K*L - 10d0*L**2
     .    - 9d0*HB*LB + 2d0*G1*LD + 6d0*G2*LD - 9d0*HB*LD - 3d0*HL*LD
     .    - 12d0*K*LD - 20d0*L*LD - 10d0*LD**2 - 3d0*HL*LL + 2d0*G1*LPP
     .    + 6d0*G2*LPP - 8d0*K*LPP - 4d0*L*LPP - 4d0*LD*LPP - 4d0*LPP**2
     .    + 4d0/3d0*G1*LT + 16d0*G3*LT - 6d0*HB*LT - 9d0*HT*LT
     .    - 3d0*LB*LT - 9d0*LD*LT - 9d0*LT**2 + 4d0/3d0*G1*LTT
     .    + 16d0*G3*LTT - 12d0*K*LTT - 6d0*L*LTT - 6d0*LD*LTT
     .    - 6d0*LTT**2 + 2d0*G1*LU + 6d0*G2*LU - 6d0*HT*LU - 8d0*K*LU
     .    - 14d0*L*LU - 6d0*LB*LU - 4d0*LD*LU - 2d0*LL*LU - 4d0*LU**2
     .    + N5*(5d0/3d0*G1**2 + 3d0*G2**2)
     .    - 18d0*DHT*DL*DLD*DLT - 15d0*DHB*DL*DLB*DLU
     .    - 5d0*DHL*DL*DLL*DLU)
     .    + DHT*DL*DLT*(4d0/3d0*G1 + 16d0*G3
     .    - 3d0*HB - 9d0*HT - 3d0*LB - 9d0*LT)))/2d0

      F(13)= (DLT*(-13d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .     + 6d0*LT + 6d0*HT + HB + LB + LD)
     .     + DL*DLD*DHT
     .     + c2*(DLT*(2743d0/162d0*G1**2 + 5d0/3d0*G1*G2 + 8d0*G2*G3
     .     + 15d0/2d0*G2**2 + 136d0/27d0*G1*G3 - 16d0/9d0*G3**2
     .     + 2d0/3d0*G1*HB - 5d0*HB**2 - HB*HL + 2d0*G1*HT + 6d0*G2*HT
     .     + 16d0*G3*HT - 5d0*HB*HT - 22d0*HT**2 - HB*L - 3d0*HT*L
     .     + 2d0/3d0*G1*LB - 10d0*HB*LB - 5d0*HT*LB - 5d0*LB**2
     .     - 4d0*HB*LD - HL*LD - 2d0*K*LD - 3d0*L*LD - 3d0*LD**2
     .     - LB*LL - 2d0*LD*LPP+ 2d0*G1*LT + 6d0*G2*LT + 16d0*G3*LT
     .     - 5d0*HB*LT - 44d0*HT*LT - 5d0*LB*LT - 3d0*LD*LT - 22d0*LT**2
     .     - 3d0*LD*LTT - 3d0*HT*LU - LB*LU - 2d0*LD*LU
     .     + N5*(65d0/27d0*G1**2 + 3d0*G2**2 + 16d0/3d0*G3**2)
     .     - 2d0*DHB*DHL*DLB*DLL - 6d0*DHT*DL*DLD*DLT
     .     - 2d0*DHB*DL*DLB*DLU)
     .     - DHT*DL*DLD*(3d0*HB + HL + 2d0*K + 3d0*L
     .     + 3d0*LD + 2d0*LPP + 3d0*LTT + 3d0*LU)
     .     - 3d0*DHB*DHT*DLB*DLD*DLU - DHL*DHT*DLD*DLL*DLU))/2d0

      F(14)= (DLB*(-7d0/9d0*G1 - 3d0*G2 - 16d0/3d0*G3
     .     + 6d0*LB + LL + LT + 6d0*HB + HT + LU)
     .     + DL*DLU*DHB + DHB*DLL*DHL
     .     + c2*(DLB*(1435d0/162d0*G1**2 + 5d0/3d0*G1*G2 + 8d0*G2*G3
     .     + 15d0/2d0*G2**2 + 40d0/27d0*G1*G3 - 16d0/9d0*G3**2
     .     + 2d0/3d0*G1*HB + 6d0*G2*HB + 16d0*G3*HB - 22d0*HB**2
     .     - 3d0*HB*HL + 4d0/3d0*G1*HT - 5d0*HB*HT - 5d0*HT**2
     .     - 3d0*HB*L - HT*L + 2d0/3d0*G1*LB + 6d0*G2*LB + 16d0*G3*LB
     .     - 44d0*HB*LB - 5d0*HT*LB - 22d0*LB**2 - 3d0*HB*LD + 2d0*G1*LL
     .     - 3d0*HL*LL - 3d0*LB*LL - 3d0*LL**2 + 4d0/3d0*G1*LT
     .     - 5d0*HB*LT - 10d0*HT*LT - 5d0*LB*LT - LD*LT - 5d0*LT**2
     .     - 4d0*HT*LU - 2d0*K*LU - 3d0*L*LU - 3d0*LB*LU - 2d0*LD*LU
     .     - 2d0*LPP*LU - 3d0*LTT*LU - 3d0*LU**2
     .     + N5*(35d0/27d0*G1**2 + 3d0*G2**2 + 16d0/3d0*G3**2)
     .     - 2d0*DHT*DL*DLD*DLT - 6d0*DHB*DL*DLB*DLU)
     .     + DHB*DHL*DLL*(2d0*G1 - 3d0*HL - 6d0*LB - 3d0*LL)
     .     - DHB*DL*DLU*(3d0*HT + 2d0*K + 3d0*L + 3d0*LD
     .     + 2d0*LPP + 3d0*LTT + 3d0*LU)
     .     - 3d0*DHB*DHT*DLD*DLT*DLU))/2d0

      F(15)= (DLL*(-3d0*G1 - 3d0*G2
     .     + 4d0*LL + 3d0*LB + 4d0*HL + LU)
     .     + DL*DLU*DHL + 3d0*DHL*DLB*DHB
     .     + c2*(DLL*(75d0/2d0*G1**2 + 3d0*G1*G2 + 15d0/2d0*G2**2
     .     + 2d0*G1*HL + 6d0*G2*HL - 9d0*HB*HL - 10d0*HL**2 - 3d0*HL*L
     .     - 2d0/3d0*G1*LB + 16d0*G3*LB - 9d0*HB*LB - 3d0*HT*LB
     .     - 9d0*LB**2 - 3d0*HL*LD + 2d0*G1*LL + 6d0*G2*LL
     .     - 20d0*HL*LL - 9d0*LB*LL - 10d0*LL**2 - 3d0*LB*LT
     .     - 3d0*HT*LU - 2d0*K*LU - 3d0*L*LU - 2d0*LD*LU - 3d0*LL*LU
     .     - 2d0*LPP*LU - 3d0*LTT*LU - 3d0*LU**2
     .     + N5*(5d0*G1**2 + 3d0*G2**2)
     .     - 18d0*DHB*DHL*DLB*DLL - 6*DHL*DL*DLL*DLU)
     .     - DHB*DHL*DLB*(2d0/3d0*G1 - 16d0*G3 + 9d0*HB + 3d0*HT
     .     + 9d0*LB + 3d0*LT)
     .     - DHL*DL*DLU*(3d0*HT + 2d0*K
     .     + 3d0*L + 3d0*LD + 2d0*LPP + 3d0*LTT)
     .     - 3d0*DHL*DHT*DLD*DLT*DLU - 3d0*LU*DLU*DHL*DL))/2d0

 111  FORMAT(5E14.6)

      END
