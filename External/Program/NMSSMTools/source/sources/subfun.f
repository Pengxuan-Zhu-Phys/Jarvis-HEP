      SUBROUTINE DIAGN(N,A,D,V,EPS)

*   Diagonalization of a real symmetric NxN matrix
*   Using the Jacobi method (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,IP,IQ,NMAX
      PARAMETER(NMAX=500)

      DOUBLE PRECISION A(N,N),D(N),V(N,N),B(N),Z(N)
      DOUBLE PRECISION EPS,SM,THR,G,H,C,S,T,THETA,TAU

      DO IP=1,N
       DO IQ=1,N
        V(IP,IQ)=0d0
       ENDDO
       V(IP,IP)=1d0
       B(IP)=A(IP,IP)
       D(IP)=B(IP)
       Z(IP)=0d0
      ENDDO

      DO I=1,NMAX

       SM=0d0
       DO IP=1,N-1
        DO IQ=IP+1,N
         SM=SM+DABS(A(IP,IQ))
        ENDDO
       ENDDO
       IF(SM.LE.EPS)RETURN

       IF(I.LT.4)THEN
        THR=.2d0*SM/N**2
       ELSE
        THR=0d0
       ENDIF

       DO IP=1,N-1
       DO IQ=IP+1,N
        IF(DABS(A(IP,IQ)).LT.EPS*MIN(DABS(D(IP)),DABS(D(IQ)))
     .     .AND. I.GT.4)THEN
         A(IP,IQ)=0d0
        ELSEIF(DABS(A(IP,IQ)).GT.THR)THEN
         H=D(IQ)-D(IP)
         IF(DABS(A(IP,IQ)).LT.EPS*DABS(H))THEN
          T=A(IP,IQ)/H
         ELSE
          THETA=.5d0*H/A(IP,IQ)
          T=1d0/(DABS(THETA)+DSQRT(1d0+THETA**2))
          IF(THETA.LT.0d0)T=-T
         ENDIF
         C=1d0/DSQRT(1d0+T**2)
         S=T*C
         TAU=S/(1d0+C)
         H=T*A(IP,IQ)
         Z(IP)=Z(IP)-H
         Z(IQ)=Z(IQ)+H
         D(IP)=D(IP)-H
         D(IQ)=D(IQ)+H
         A(IP,IQ)=0d0
         DO J=1,IP-1
          G=A(J,IP)
          H=A(J,IQ)
          A(J,IP)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
         ENDDO
         DO J=IP+1,IQ-1
          G=A(IP,J)
          H=A(J,IQ)
          A(IP,J)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
         ENDDO
         DO J=IQ+1,N
          G=A(IP,J)
          H=A(IQ,J)
          A(IP,J)=G-S*(H+G*TAU)
          A(IQ,J)=H+S*(G-H*TAU)
         ENDDO
         DO J=1,N
          G=V(J,IP)
          H=V(J,IQ)
          V(J,IP)=G-S*(H+G*TAU)
          V(J,IQ)=H+S*(G-H*TAU)
         ENDDO
        ENDIF
       ENDDO
       ENDDO

       DO IP=1,N
        B(IP)=B(IP)+Z(IP)
        D(IP)=B(IP)
        Z(IP)=0d0
       ENDDO

      ENDDO

      RETURN
      END


      SUBROUTINE SORTN(N,D,V)

*   Reordering of the eigenvalues D(I), I=1..N
*   and corresponding eigenvectors V(J,I), J=1..N in ascending order
*   (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,K
      DOUBLE PRECISION D(N),V(N,N),P

      DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1,N
        IF(D(J).LT.P)THEN
         K=J
         P=D(J)
        ENDIF
       ENDDO
       IF(K.NE.I)THEN
        D(K)=D(I)
        D(I)=P
        DO J=1,N
         P=V(J,I)
         V(J,I)=V(J,K)
         V(J,K)=P
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END


      SUBROUTINE SORTNA(N,D,V)

*   Reordering of the absolute value of the eigenvalues D(I), I=1..N
*   and corresponding eigenvectors V(J,I), J=1..N in ascending order
*   (from Numerical Recipes)

      IMPLICIT NONE

      INTEGER N,I,J,K
      DOUBLE PRECISION D(N),V(N,N),P

      DO I=1,N-1
       K=I
       P=D(I)
       DO J=I+1,N
        IF(DABS(D(J)).LT.DABS(P))THEN
         K=J
         P=D(J)
        ENDIF
       ENDDO
       IF(K.NE.I)THEN
        D(K)=D(I)
        D(I)=P
        DO J=1,N
         P=V(J,I)
         V(J,I)=V(J,K)
         V(J,K)=P
        ENDDO
       ENDIF
      ENDDO

      RETURN
      END


      DOUBLE PRECISION FUNCTION RUNM(Q,NF)

*   Subroutine to calculate the quark running masses

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      PARAMETER (NN=6)
      PARAMETER (ZETA3 = 1.202056903159594d0)

      DIMENSION AM(NN),YMSB(NN)

      COMMON/ALS/XLAMBDA,AMCA,AMBA,AMTA,N0A
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
      COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

      B0(NF)= (33d0-2d0*NF)/12d0
      B1(NF)= (102d0-38d0/3d0*NF)/16d0
      B2(NF)= (2857d0/2d0-5033d0/18d0*NF+325d0/54d0*NF**2)/64d0
      G0(NF)= 1d0
      G1(NF)= (202d0/3d0-20d0/9d0*NF)/16d0
      G2(NF)= (1249d0-(2216d0/27d0+160d0/3d0*ZETA3)*NF
     .      - 140d0/81d0*NF**2)/64d0
      C1(NF)= G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2
      C2(NF)= ((G1(NF)/B0(NF) - B1(NF)*G0(NF)/B0(NF)**2)**2
     .      + G2(NF)/B0(NF) + B1(NF)**2*G0(NF)/B0(NF)**3
     .      - B1(NF)*G1(NF)/B0(NF)**2 - B2(NF)*G0(NF)/B0(NF)**2)/2d0
      TRAN(X,XK)= 1d0+4d0/3d0*ALPHAS(X,2)/PI+XK*(ALPHAS(X,2)/PI)**2
      CQ(X,NF)= (2d0*B0(NF)*X)**(G0(NF)/B0(NF))
     .  * (1d0+C1(NF)*X+C2(NF)*X**2)

      PI= 4d0*DATAN(1d0)
      ACC= 1d-8
      AM(1)= 0
      AM(2)= 0
      AM(3)= AMS
      AM(4)= AMC
      AM(5)= AMBP
      AM(6)= AMT
      XK= 16.11d0
      DO 1 I=1,NF-1
       XK= XK - 1.04d0*(1d0-AM(I)/AM(NF))
1     CONTINUE
      IF(NF.GE.4)THEN
       XMSB= AM(NF)/TRAN(AM(NF),0d0)
       XMHAT= XMSB/CQ(ALPHAS(AM(NF),2)/PI,NF)
      ELSE
       XMSB= 0
       XMHAT= 0
      ENDIF
      YMSB(3)= AMS
      IF(NF.EQ.3)THEN
       YMSB(4)= YMSB(3)*CQ(ALPHAS(AM(4),2)/PI,3)/
     .            CQ(ALPHAS(1d0,2)/PI,3)
       YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .            CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.4)THEN
       YMSB(4)= XMSB
       YMSB(5)= YMSB(4)*CQ(ALPHAS(AM(5),2)/PI,4)/
     .            CQ(ALPHAS(AM(4),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.5)THEN
       YMSB(5)= XMSB
       YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .            CQ(ALPHAS(AM(5),2)/PI,4)
       YMSB(6)= YMSB(5)*CQ(ALPHAS(AM(6),2)/PI,5)/
     .            CQ(ALPHAS(AM(5),2)/PI,5)
      ELSEIF(NF.EQ.6)THEN
       YMSB(6)= XMSB
       YMSB(5)= YMSB(6)*CQ(ALPHAS(AM(5),2)/PI,5)/
     .            CQ(ALPHAS(AM(6),2)/PI,5)
       YMSB(4)= YMSB(5)*CQ(ALPHAS(AM(4),2)/PI,4)/
     .            CQ(ALPHAS(AM(5),2)/PI,4)
      ENDIF
      IF(Q.LT.AMC)THEN
       N0= 3
       Q0= 1d0
      ELSEIF(Q.LE.AMBP)THEN
       N0= 4
       Q0= AMC
      ELSEIF(Q.LE.AMT)THEN
       N0= 5
       Q0= AMBP
      ELSE
       N0= 6
       Q0= AMT
      ENDIF
      IF(NF.GT.3)THEN
       XKFAC= TRAN(AM(NF),0d0)/TRAN(AM(NF),XK)
      ELSE
       XKFAC= 1d0
      ENDIF
      RUNM= YMSB(N0)*CQ(ALPHAS(Q,2)/PI,N0)/
     .         CQ(ALPHAS(Q0,2)/PI,N0)
     .       * XKFAC

      RETURN
      END


*  Running alpha_s and aux. subroutines/functions as in HDECAY

      DOUBLE PRECISION FUNCTION ALPHAS(Q,N)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      DIMENSION XLB(6)

      COMMON/ALSLAM/XLB1(6),XLB2(6),XLB3(6)
      COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      B2(NF)=27d0/2d0*(2857d0-5033d0/9d0*NF+325d0/27d0*NF**2)/B0(NF)**3
      ALS1(NF,X)=12d0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
      ALS2(NF,X)=12d0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .    *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .     /DLOG(X**2/XLB(NF)**2))
      ALS3(NF,X)=12d0*PI/(B0(NF)*DLOG(X**2/XLB(NF)**2))
     .          *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB(NF)**2))
     .           /DLOG(X**2/XLB(NF)**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB(NF)**2))**2
     .                      -DLOG(DLOG(X**2/XLB(NF)**2))-1d0)+B2(NF))
     .           /DLOG(X**2/XLB(NF)**2)**2)

      PI=4d0*DATAN(1d0)

      IF(N.EQ.1)THEN
       DO 1 I=1,6
        XLB(I)=XLB1(I)
1     CONTINUE
      ELSE
       DO 2 I=1,6
        XLB(I)=XLB2(I)
2     CONTINUE
      ENDIF
      IF(Q.LT.AMC)THEN
       NF=3
      ELSEIF(Q.LE.AMBP)THEN
       NF=4
      ELSEIF(Q.LE.AMT)THEN
       NF=5
      ELSE
       NF=6
      ENDIF
      IF(N.EQ.1)THEN
        ALPHAS=ALS1(NF,Q)
      ELSEIF(N.EQ.2)THEN
        ALPHAS=ALS2(NF,Q)
      ELSE
        ALPHAS=ALS3(NF,Q)
      ENDIF

      RETURN
      END


      SUBROUTINE ALSINI(ACC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      DIMENSION XLB(6)

      COMMON/ALSLAM/XLB1(6),XLB2(6),XLB3(6)
      COMMON/ALS/XLAMBDA,AMC,AMBP,AMT,N0

      PI=4d0*DATAN(1d0)
      XLB1(1)=0d0
      XLB1(2)=0d0
      XLB2(1)=0d0
      XLB2(2)=0d0
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2d0/25d0)
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2d0/23d0)
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
      ENDIF
      DO 1 I=1,6
       XLB1(I)=XLB(I)
1     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2d0/25d0)
     .       *(2d0*DLOG(AMC/XLB(3)))**(-107d0/1875d0)
       XLB(4)=XITER(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .       *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .       *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .       *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .      *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2d0/23d0)
     .      *(2d0*DLOG(AMT/XLB(6)))**(321d0/3703d0)
       XLB(5)=XITER(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .       *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .       *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 2 I=1,6
       XLB2(I)=XLB(I)
2     CONTINUE
      IF(N0.EQ.3)THEN
       XLB(3)=XLAMBDA
       XLB(4)=XLB(3)*(XLB(3)/AMC)**(2d0/25d0)
     .             *(2d0*DLOG(AMC/XLB(3)))**(-107d0/1875d0)
       XLB(4)=XITER3(AMC,XLB(3),3,XLB(4),4,ACC)
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .             *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER3(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .            *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.4)THEN
       XLB(4)=XLAMBDA
       XLB(5)=XLB(4)*(XLB(4)/AMBP)**(2d0/23d0)
     .             *(2d0*DLOG(AMBP/XLB(4)))**(-963d0/13225d0)
       XLB(5)=XITER3(AMBP,XLB(4),4,XLB(5),5,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .             *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .            *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.5)THEN
       XLB(5)=XLAMBDA
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .             *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER3(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .             *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
       XLB(6)=XLB(5)*(XLB(5)/AMT)**(2d0/21d0)
     .            *(2d0*DLOG(AMT/XLB(5)))**(-321d0/3381d0)
       XLB(6)=XITER3(AMT,XLB(5),5,XLB(6),6,ACC)
      ELSEIF(N0.EQ.6)THEN
       XLB(6)=XLAMBDA
       XLB(5)=XLB(6)*(XLB(6)/AMT)**(-2d0/23d0)
     .            *(2d0*DLOG(AMT/XLB(6)))**(321d0/3703d0)
       XLB(5)=XITER3(AMT,XLB(6),6,XLB(5),5,ACC)
       XLB(4)=XLB(5)*(XLB(5)/AMBP)**(-2d0/25d0)
     .             *(2d0*DLOG(AMBP/XLB(5)))**(963d0/14375d0)
       XLB(4)=XITER3(AMBP,XLB(5),5,XLB(4),4,ACC)
       XLB(3)=XLB(4)*(XLB(4)/AMC)**(-2d0/27d0)
     .             *(2d0*DLOG(AMC/XLB(4)))**(107d0/2025d0)
       XLB(3)=XITER3(AMC,XLB(4),4,XLB(3),3,ACC)
      ENDIF
      DO 3 I=3,6
       XLB3(I)=XLB(I)
3     CONTINUE

      RETURN
      END


      DOUBLE PRECISION FUNCTION XITER(Q,XLB1,NF1,XLB,NF2,ACC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      ALS2(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .        *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .        /DLOG(X**2/XLB**2))
      AA(NF)=12d0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      XIT(A,B,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)))
      PI=4d0*DATAN(1d0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      ALP=ALS2(NF1,Q,XLB1)
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      XX=XIT(A,B,X)
      XLB2=Q*DEXP(-XX/2d0)
      Y1=ALS2(NF1,Q,XLB1)
      Y2=ALS2(NF2,Q,XLB2)
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
      XITER=XLB2

      RETURN
      END


      DOUBLE PRECISION FUNCTION XITER3(Q,XLB1,NF1,XLB,NF2,ACC)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      B2(NF)=27d0/2d0*(2857d0-5033d0/9d0*NF+325d0/27d0*NF**2)/B0(NF)**3
      ALS3(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1d0)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2)
      AA(NF)=12d0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      XIT(A,B,C,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)
     .          *(1d0-(A*B*(DLOG(X)**2-DLOG(X)-1d0)+C/B)/X/DLOG(X))))
      PI=4d0*DATAN(1d0)
      XLB2=XLB
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB2**2)
      IF(NF1.LT.NF2)THEN
       DELTA = 7d0*ALS3(NF1,Q,XLB1)**2/PI**2/24d0
       ALP=ALS3(NF1,Q,XLB1)*(1d0+DELTA)
      ELSE
       DELTA = 7d0*ALS3(NF1,Q,XLB1)**2/PI**2/24d0
       ALP=ALS3(NF1,Q,XLB1)/(1d0+DELTA)
      ENDIF
      A=AA(NF2)/ALP
      B=BB(NF2)*ALP
      C=CC(NF2)*ALP
      XX=XIT(A,B,C,X)
      XLB2=Q*DEXP(-XX/2d0)
      IF(NF1.LT.NF2)THEN
       DELTA = 7d0*ALS3(NF1,Q,XLB1)**2/PI**2/24d0
       Y1=ALS3(NF1,Q,XLB1)*(1d0+DELTA)
       Y2=ALS3(NF2,Q,XLB2)
      ELSE
       DELTA = 7d0*ALS3(NF1,Q,XLB1)**2/PI**2/24d0
       Y1=ALS3(NF1,Q,XLB1)/(1d0+DELTA)
       Y2=ALS3(NF2,Q,XLB2)
      ENDIF
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
       XITER3=XLB2

      RETURN
      END


      DOUBLE PRECISION FUNCTION XITLA(NO,ALP,ACC)

*  Iteration routine to determine improved Lambda's

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,gg1,gg2,S2TW
      COMMON/SMSPEC/AMS,AMC,AMB,AMBP,AMT,AMTAU,AMMUON,AMZ,AMW

      B0(NF)=33d0-2d0*NF
      B1(NF)=6d0*(153d0-19d0*NF)/B0(NF)**2
      B2(NF)=27d0/2d0*(2857d0-5033d0/9d0*NF+325d0/27d0*NF**2)/B0(NF)**3
      ALS2(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .        *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .        /DLOG(X**2/XLB**2))
      ALS3(NF,X,XLB)=12d0*PI/(B0(NF)*DLOG(X**2/XLB**2))
     .          *(1d0-B1(NF)*DLOG(DLOG(X**2/XLB**2))
     .           /DLOG(X**2/XLB**2)
     .           +(B1(NF)**2*(DLOG(DLOG(X**2/XLB**2))**2
     .                      -DLOG(DLOG(X**2/XLB**2))-1d0)+B2(NF))
     .           /DLOG(X**2/XLB**2)**2)
      AA(NF)=12d0*PI/B0(NF)
      BB(NF)=B1(NF)/AA(NF)
      CC(NF)=B2(NF)/AA(NF)
      XIT(A,B,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)))
      XIT3(A,B,C,X)=A/2d0*(1d0+DSQRT(1d0-4d0*B*DLOG(X)
     .          *(1d0-(A*B*(DLOG(X)**2-DLOG(X)-1d0)+C/B)/X/DLOG(X))))
      PI=4d0*DATAN(1d0)
      NF=5
      Q=AMZ
      XLB=Q*DEXP(-AA(NF)/ALP/2d0)
      IF(NO.EQ.1)GOTO 111
      II=0
1     II=II+1
      X=DLOG(Q**2/XLB**2)
      A=AA(NF)/ALP
      B=BB(NF)*ALP
      C=CC(NF)*ALP
      IF(NO.EQ.2)THEN
       XX=XIT(A,B,X)
      ELSE
       XX=XIT3(A,B,C,X)
      ENDIF
      XLB=Q*DEXP(-XX/2d0)
      Y1=ALP
      IF(NO.EQ.2)THEN
       Y2=ALS2(NF,Q,XLB)
      ELSE
       Y2=ALS3(NF,Q,XLB)
      ENDIF
      DY=DABS(Y2-Y1)/Y1
      IF(DY.GE.ACC) GOTO 1
111   XITLA=XLB

      RETURN
      END


      DOUBLE PRECISION FUNCTION FINT(Z,XX,YY)

*  One-dimensional cubic interpolation
*  Z  = wanted point
*  XX = array of 4 discrete x-values around Z
*  YY = array of 4 discrete function-values around Z

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION XX(4),YY(4)

      X = DLOG(Z)
      X0=DLOG(XX(1))
      X1=DLOG(XX(2))
      X2=DLOG(XX(3))
      X3=DLOG(XX(4))
      Y0=DLOG(YY(1))
      Y1=DLOG(YY(2))
      Y2=DLOG(YY(3))
      Y3=DLOG(YY(4))
      A0=(X-X1)*(X-X2)*(X-X3)/(X0-X1)/(X0-X2)/(X0-X3)
      A1=(X-X0)*(X-X2)*(X-X3)/(X1-X0)/(X1-X2)/(X1-X3)
      A2=(X-X0)*(X-X1)*(X-X3)/(X2-X0)/(X2-X1)/(X2-X3)
      A3=(X-X0)*(X-X1)*(X-X2)/(X3-X0)/(X3-X1)/(X3-X2)
      FINT=DEXP(A0*Y0+A1*Y1+A2*Y2+A3*Y3)

      RETURN
      END


*   Spence function and auxiliary functions as in HDECAY

      DOUBLE PRECISION FUNCTION SP(X)

*  REAL dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE COMPLEX CX,LI2

      CX = DCMPLX(X,0d0)
      SP = DREAL(LI2(CX))

      RETURN
      END


      DOUBLE COMPLEX FUNCTION LI2(X)

*  COMPLEX dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE COMPLEX X,Y,CLI2

      COMMON/CONST/ZETA2,ZETA3

      ZERO=1d-16
      XR=DREAL(X)
      XI=DIMAG(X)
      R2=XR*XR+XI*XI
      LI2=0
      IF(R2.LE.ZERO)THEN
        LI2=X
        RETURN
      ENDIF
      RR=XR/R2
      IF(R2.EQ.1d0.AND.XI.EQ.0d0)THEN
        IF(XR.EQ.1d0)THEN
          LI2=DCMPLX(ZETA2)
        ELSE
          LI2=-DCMPLX(ZETA2/2d0)
        ENDIF
        RETURN
      ELSEIF(R2.GT.1d0.AND.RR.GT.0.5d0)THEN
        Y=(X-1d0)/X
        LI2=CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1d0-X)+0.5d0*CDLOG(X)**2
        RETURN
      ELSEIF(R2.GT.1d0.AND.RR.LE.0.5d0)THEN
        Y=1d0/X
        LI2=-CLI2(Y)-ZETA2-0.5d0*CDLOG(-X)**2
        RETURN
      ELSEIF(R2.LE.1d0.AND.XR.GT.0.5d0)THEN
        Y=1d0-X
        LI2=-CLI2(Y)+ZETA2-CDLOG(X)*CDLOG(1d0-X)
       RETURN
      ELSEIF(R2.LE.1d0.AND.XR.LE.0.5d0)THEN
        Y=X
        LI2=CLI2(Y)
        RETURN

      ENDIF
      END


      DOUBLE COMPLEX FUNCTION CLI2(X)

*  Taylor-expansion for complex dilogarithm (Spence-function)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DOUBLE COMPLEX X,Z

      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/POLY/NBER

      N=NBER-1
      Z=-CDLOG(1d0-X)
      CLI2=B2(NBER)
      DO 111 I=N,1,-1
        CLI2=Z*CLI2+B2(I)
111   CONTINUE
      CLI2=Z**2*CLI2+Z

      RETURN
      END


      DOUBLE PRECISION FUNCTION FACULT(N)

*  DOUBLE PRECISION version of FACULTY

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

      FACULT=1d0
      IF(N.EQ.0)RETURN
      DO 999 I=1,N
        FACULT=FACULT*DFLOAT(I)
999   CONTINUE

      RETURN
      END


      SUBROUTINE BERNINI(N)

*  Initialization of coefficients for polylogarithms

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION B(18),PB(19)

      COMMON/BERNOULLI/B2(18),B12(18),B3(18)
      COMMON/CONST/ZETA2,ZETA3
      COMMON/POLY/NBER

      NBER=N
      PI=4d0*DATAN(1d0)

      B(1)=-1d0/2d0
      B(2)=1d0/6d0
      B(3)=0d0
      B(4)=-1d0/30d0
      B(5)=0d0
      B(6)=1d0/42d0
      B(7)=0d0
      B(8)=-1d0/30d0
      B(9)=0d0
      B(10)=5d0/66d0
      B(11)=0d0
      B(12)=-691d0/2730d0
      B(13)=0d0
      B(14)=7d0/6d0
      B(15)=0d0
      B(16)=-3617d0/510d0
      B(17)=0d0
      B(18)=43867d0/798d0
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0

      DO 995 I=1,18
        B2(I)=B(I)/FACULT(I+1)
        B12(I)=DFLOAT(I+1)/FACULT(I+2)*B(I)/2d0
        PB(I+1)=B(I)
        B3(I)=0d0
995   CONTINUE
      PB(1)=1d0
      DO 996 I=1,18
      DO 996 J=0,I
       B3(I)=B3(I)+PB(J+1)*PB(I-J+1)/FACULT(I-J)/FACULT(J+1)
     .       /DFLOAT(I+1)
996   CONTINUE

      RETURN
      END


*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*   Passarino-Veltman one- and two-points functions A0, B0 and B1
*   orig from LoopTools, http://www.feynarts.de/looptools/
*   taken from Suspect2.3, modified by S. Kraml, 7 March 2005
*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      double precision function NMA0(m2,q)

      implicit none
      double precision m2,q
      if(m2.ne.0d0) then
       NMA0 = m2 * (1d0-dlog( m2/q ))
      else
       NMA0 = 0d0
      endif
      end

      double precision function NMB0(p,m1,m2,q)

*     note: all input is quadratical: p=p^2, m1=m1^2, m2=m2^2, q=q^2

      implicit none
      double precision p, m1, m2
      double precision mudim2, divergence, lambda2, q
      double precision acc, eps, minacc
      double complex x1, x2, y1, y2, r, be0
      double complex Ieps, onePeps, oneMeps
      COMMON/cutoff/mudim2, divergence, lambda2
      parameter (acc = 1d-12)
      parameter (eps = 1d-20)
      parameter (Ieps = (0d0,1d0)*eps)
      parameter (onePeps = 1d0 + Ieps)
      parameter (oneMeps = 1d0 - Ieps)

      double complex fpv, xlogx
      external fpv, xlogx

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

      NMB0 = dble(be0 + divergence)

      end


      double precision function NMB1(s,mi,mj,q)

* note: all input is quadratical: s=p^2, mi=m1^2, mj=m2^2, q=q^2

      implicit none
      double precision s,mi,mj,NMB0,NMA0,q

      if(mi.eq.mj) then
       NMB1 = NMB0(s,mi,mj,q)/2d0
      else
       NMB1= (NMA0(mj,q) - NMA0(mi,q)
     .   + (s+mi-mj)*NMB0(s,mi,mj,q))/(2d0*s)
      endif
      end


*---------------------------------------------------------------------
* auxiliary functions used by the B0,B1 two-point functions
* from Looptools http://www.feynarts.de/looptools/
*---------------------------------------------------------------------

      subroutine roots(p, m1, m2, x1, x2, y1, y2, r)

      implicit none
      double precision p, m1, m2
      double complex x1, x2, y1, y2, r
      double precision mudim2, divergence, lambda2
      COMMON/cutoff/mudim2, divergence, lambda2
      double precision acc, eps
      double complex Ieps, onePeps, oneMeps
      parameter (acc = 1d-12)
      parameter (eps = 1d-20)
      parameter (Ieps = (0d0,1d0)*eps)
      parameter (onePeps = 1d0 + Ieps)
      parameter (oneMeps = 1d0 - Ieps)
      double precision q

      r = sqrt(dcmplx(p*(p - 2*(m1 + m2)) + (m1 - m2)**2))
      q = p + m1 - m2
      x1 = (q + r)/2d0/p
      x2 = (q - r)/2d0/p
      if(abs(x2) .gt. abs(x1)) then
        x1 = m1/p/x2
      else if(abs(x1) .gt. abs(x2)) then
        x2 = m1/p/x1
      endif
      x1 = x1 + abs(p*x1)/p*Ieps
      x2 = x2 - abs(p*x2)/p*Ieps
      q = p - m1 + m2
      y2 = (q + r)/2d0/p
      y1 = (q - r)/2d0/p
      if(abs(y2) .gt. abs(y1)) then
        y1 = m2/p/y2
      else if(abs(y1) .gt. abs(y2)) then
        y2 = m2/p/y1
      endif
      y1 = y1 - abs(p*y1)/p*Ieps
      y2 = y2 + abs(p*y2)/p*Ieps
      end


      double complex function fpv(n, x, y)

      implicit none
      integer n
      double complex x, y
      double precision mudim2, divergence, lambda2
      COMMON/cutoff/mudim2, divergence, lambda2
      double precision acc, eps
      double complex Ieps, onePeps, oneMeps
      parameter (acc = 1d-12)
      parameter (eps = 1d-20)
      parameter (Ieps = (0,1)*eps)
      parameter (onePeps = 1 + Ieps)
      parameter (oneMeps = 1 - Ieps)
      integer m
      double complex xm
      if(abs(x) .lt. 10d0) then
        if(n .eq. 0) then
          fpv = -log(-y/x)
        else if(abs(x) .lt. acc) then
          fpv = -1d0/n
        else
          fpv = 0
          xm = 1
          do m = 0, n - 1
            fpv = fpv - xm/(n - m)
            xm = xm*x
          enddo
          fpv = fpv - xm*log(-y/x)
        endif
      else
        fpv = 0
        xm = 1
        do m = 1, 30
          xm = xm/x
          fpv = fpv + xm/(m + n)
          if(abs(xm/fpv) .lt. acc**2) return
        enddo
      endif
      end


      double complex function yfpv(n, x, y)

      implicit none
      integer n
      double complex x, y
      double complex fpv
      external fpv
      if(abs(y) .eq. 0d0) then
        yfpv = 0
      else
        yfpv = y*fpv(n, x, y)
      endif
      end


      double complex function xlogx(x)

      implicit none
      double complex x
      if(abs(x) .eq. 0d0) then
        xlogx = 0
      else
        xlogx = x*log(x)
      endif
      end


      DOUBLE PRECISION FUNCTION RUNMB(Q)
      
*   Subroutine to calculate the running b quark mass for Q > MB

      IMPLICIT NONE
      DOUBLE PRECISION Q
      DOUBLE PRECISION PI,ALPHAS,ALMB,ALMT,ALQ,U5MTMB,U6QMT
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW

      PI=4d0*DATAN(1d0)

      ALQ=ALPHAS(Q,2)
      ALMB=ALPHAS(MB,2)

      IF(Q.LE.MT) THEN

       RUNMB=MB*(ALQ/ALMB)**(12d0/23d0)*(1d0+7462d0*(ALQ-ALMB)/
     .      (4d0*PI*1587d0))

      ELSE

       ALMT=ALPHAS(MT,2)
       U5MTMB=(ALMT/ALMB)**(12d0/23d0)*(1d0+7462d0*(ALMT-ALMB)/
     .      (4d0*PI*1587d0))
        U6QMT=(ALQ/ALMT)**(4d0/7d0)*(1d0+7398d0*(ALQ-ALMT)/
     .       (4d0*PI*1323d0))
       RUNMB=MB*U6QMT*U5MTMB

      ENDIF

      END


      DOUBLE PRECISION FUNCTION RAN2(IDUM)

      IMPLICIT NONE

      INTEGER IDUM
      DOUBLE PRECISION DA,DB,DC
      PARAMETER(DA=16807d0,DB=2147483647d0,DC=2147483648d0)

      IDUM=INT(ABS(MOD(DA*IDUM,DB)+0.5d0))
      RAN2=DFLOAT(IDUM)/DC

      END


      DOUBLE PRECISION FUNCTION GAU(IDUM)

      IMPLICIT NONE

      INTEGER IDUM,ISET
      DOUBLE PRECISION F,SET,R,V1,V2,RAN2
      SAVE ISET,SET
      DATA ISET/0/
      IF(ISET.EQ.0)THEN
 5      V1=2d0*RAN2(IDUM)-1d0
        V2=2d0*RAN2(IDUM)-1d0
        R=V1**2+V2**2
        IF(R.GE.1d0.OR.R.EQ.0d0)GOTO 5
        F=DSQRT(-2d0*DLOG(R)/R)
        SET=V1*F
        GAU=V2*F
        ISET=1
      ELSE
        GAU=SET
        ISET=0
      ENDIF
      END


      DOUBLE PRECISION FUNCTION UPSILONTAU(MX)

*  BABAR constraints on BR(Y -> A gamma)*BR(A -> tau tau) from 1210.5669 Fig. 3

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=127)
      DOUBLE PRECISION MX,X(N),M(N)

      DATA M/3.5964d0,3.7003d0,3.7993d0,3.8983d0,3.9991d0,4.0967d0,
     .4.1629d0,4.2377d0,4.3060d0,4.3669d0,4.4333d0,4.4912d0,4.5557d0,
     .4.6196d0,4.6775d0,4.7417d0,4.7865d0,4.8927d0,4.9544d0,5.0053d0,
     .5.0581d0,5.1134d0,5.1598d0,5.2079d0,5.2690d0,5.3038d0,5.3585d0,
     .5.4450d0,5.4903d0,5.6329d0,5.7524d0,5.7960d0,5.8934d0,5.9156d0,
     .5.9570d0,5.9898d0,6.0223d0,6.1035d0,6.1730d0,6.2408d0,6.2867d0,
     .6.3131d0,6.3576d0,6.3878d0,6.4259d0,6.5166d0,6.5608d0,6.5900d0,
     .6.6161d0,6.6514d0,6.6778d0,6.7112d0,6.7707d0,6.7929d0,6.8260d0,
     .6.8569d0,6.9183d0,6.9447d0,6.9616d0,6.9989d0,7.0514d0,7.0753d0,
     .7.1075d0,7.1336d0,7.1569d0,7.1853d0,7.2050d0,7.2423d0,7.2931d0,
     .7.3599d0,7.3885d0,7.4082d0,7.4327d0,7.4877d0,7.5667d0,7.5945d0,
     .7.6142d0,7.6323d0,7.7238d0,7.7438d0,7.7638d0,7.8900d0,7.9101d0,
     .7.9270d0,7.9690d0,7.9884d0,8.0040d0,8.0240d0,8.0412d0,8.0584d0,
     .8.1313d0,8.1516d0,8.1688d0,8.1819d0,8.2197d0,8.2308d0,8.2553d0,
     .8.2839d0,8.3103d0,8.3275d0,8.3842d0,8.4042d0,8.4221d0,8.4357d0,
     .8.4557d0,8.5091d0,8.5199d0,8.5455d0,8.6712d0,8.6848d0,8.7090d0,
     .8.7173d0,8.7513d0,8.7757d0,8.8133d0,8.8266d0,8.9067d0,8.9201d0,
     .8.9392d0,8.9851d0,9.0404d0,9.0643d0,9.0815d0,9.0951d0,9.1073d0,
     .9.1157d0,9.1935d0/

      DATA X/9.6363d-6,1.1592d-5,1.6331d-5,1.9071d-5,1.1368d-5,
     .8.9858d-6,1.0457d-5,1.0633d-5,9.7667d-6,9.1600d-6,9.8759d-6,
     .1.0707d-5,1.1254d-5,1.1501d-5,1.1062d-5,9.4220d-6,8.9810d-6,
     .8.9456d-6,8.9933d-6,9.2746d-6,9.3759d-6,9.2723d-6,9.2723d-6,
     .9.1065d-6,9.4130d-6,9.5513d-6,9.6692d-6,9.5092d-6,9.7860d-6,
     .1.1835d-5,1.3486d-5,1.4591d-5,1.8993d-5,2.0172d-5,2.1389d-5,
     .2.2304d-5,2.3130d-5,2.3755d-5,2.5095d-5,2.5981d-5,2.7006d-5,
     .2.7555d-5,2.6869d-5,2.4955d-5,2.0875d-5,1.1637d-5,9.6551d-6,
     .9.0927d-6,8.7316d-6,8.8129d-6,8.9733d-6,1.0202d-5,1.4303d-5,
     .1.5713d-5,1.6727d-5,1.6727d-5,1.5698d-5,1.5785d-5,1.6017d-5,
     .1.6107d-5,1.5345d-5,1.5116d-5,1.5502d-5,1.6429d-5,1.6970d-5,
     .1.7501d-5,1.7501d-5,1.6900d-5,1.5365d-5,1.3983d-5,1.3695d-5,
     .1.3695d-5,1.4626d-5,1.9321d-5,3.5421d-5,3.7549d-5,3.6308d-5,
     .3.2669d-5,9.6585d-6,9.3982d-6,9.8574d-6,1.6491d-5,1.7682d-5,
     .1.8783d-5,2.4249d-5,2.5646d-5,3.7287d-5,4.2060d-5,4.3951d-5,
     .4.3305d-5,1.8772d-5,1.6894d-5,1.6742d-5,1.7322d-5,2.0555d-5,
     .2.1070d-5,2.1458d-5,2.1031d-5,2.1256d-5,2.2160d-5,3.9401d-5,
     .4.0850d-5,3.7669d-5,3.6663d-5,4.1364d-5,7.8771d-5,8.2049d-5,
     .7.9985d-5,2.8887d-5,2.7497d-5,2.7497d-5,2.8245d-5,4.8374d-5,
     .4.7661d-5,1.9951d-5,2.1233d-5,1.0924d-4,1.1849d-4,1.2020d-4,
     .5.3685d-5,4.0960d-5,3.0863d-5,3.0255d-5,3.1100d-5,3.1281d-5,
     .3.0213d-5,1.2731d-4/

      UPSILONTAU=0d0

      IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

      DO I=2,N
       IF(MX.LT.M(I))THEN
        UPSILONTAU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION UPSILONMU(MX)

*  BABAR constraints on BR(Y -> A gamma)*BR(A -> mu mu) from 1210.0287, Fig. 5(c)

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=208)
      DOUBLE PRECISION MX,X(N),M(N)

      DATA M/1.9811d-1,2.2708d-1,2.4048d-1,2.5986d-1,2.8142d-1,
     .2.9356d-1,3.1166d-1,3.2506d-1,3.4098d-1,3.6146d-1,3.7831d-1,
     .3.9769d-1,4.1109d-1,4.2683d-1,4.4621d-1,4.6526d-1,4.8226d-1,
     .4.9692d-1,5.1486d-1,5.3315d-1,5.5126d-1,5.6829d-1,5.8024d-1,
     .6.0089d-1,6.1901d-1,6.3601d-1,6.5430d-1,6.6753d-1,6.8453d-1,
     .7.0139d-1,7.2439d-1,7.3887d-1,7.5336d-1,7.7510d-1,7.9322d-1,
     .8.0770d-1,8.2341d-1,8.4279d-1,8.5729d-1,8.7430d-1,8.9133d-1,
     .9.0944d-1,9.2501d-1,9.4422d-1,9.6706d-1,9.8138d-1,9.9570d-1,
     .1.0125d0,1.0176d0,1.0310d0,1.0455d0,1.0625d0,1.0857d0,1.0965d0,
     .1.1174d0,1.1355d0,1.1498d0,1.1655d0,1.1871d0,1.1991d0,1.2268d0,
     .1.2376d0,1.2570d0,1.2740d0,1.2889d0,1.3034d0,1.3251d0,1.3385d0,
     .1.3582d0,1.3714d0,1.3919d0,1.4077d0,1.4270d0,1.4441d0,1.4607d0,
     .1.4790d0,1.4960d0,1.5080d0,1.5301d0,1.5471d0,1.5628d0,1.5759d0,
     .1.5917d0,1.6123d0,1.6342d0,1.6485d0,1.6835d0,1.6944d0,1.7186d0,
     .1.7353d0,1.7570d0,1.7677d0,1.7992d0,1.8221d0,1.8366d0,1.8655d0,
     .1.8802d0,1.8972d0,1.9166d0,1.9408d0,1.9590d0,1.9709d0,1.9941d0,
     .2.0086d0,2.0265d0,2.0435d0,2.0617d0,2.0727d0,2.0955d0,2.1174d0,
     .2.1281d0,2.1464d0,2.1618d0,2.1812d0,2.1995d0,2.2152d0,2.2333d0,
     .2.2513d0,2.2683d0,2.2806d0,2.3147d0,2.3375d0,2.3543d0,2.3665d0,
     .2.3869d0,2.4088d0,2.4208d0,2.4365d0,2.4606d0,2.4701d0,2.4918d0,
     .2.5038d0,2.5279d0,2.5375d0,2.5666d0,2.5775d0,2.5920d0,2.6088d0,
     .2.6293d0,2.6476d0,2.6632d0,2.6764d0,2.6967d0,2.7183d0,2.7315d0,
     .2.7509d0,2.7630d0,2.7824d0,2.7992d0,2.8173d0,2.8355d0,2.8525d0,
     .2.8693d0,2.8825d0,2.9044d0,2.9200d0,2.9380d0,2.9550d0,2.9731d0,
     .2.9887d0,3.0044d0,3.0238d0,3.0359d0,3.0591d0,3.0917d0,3.1075d0,
     .3.1256d0,3.1461d0,3.1776d0,3.1908d0,3.2102d0,3.2259d0,3.2609d0,
     .3.2803d0,3.3141d0,3.3322d0,3.3480d0,3.3648d0,3.3828d0,3.3996d0,
     .3.4178d0,3.4335d0,3.4563d0,3.4679d0,3.4851d0,3.5336d0,3.5593d0,
     .3.5908d0,3.6247d0,3.6634d0,3.6779d0,3.6900d0,3.7082d0,3.7263d0,
     .3.7408d0,3.7757d0,3.8050d0,3.8148d0,3.8280d0,3.8499d0,3.8632d0,
     .3.8814d0,3.9008d0,3.9189d0,3.9334d0,3.9492d0,3.9662d0,3.9941d0/

      DATA X/2.8268d-7,2.8425d-7,3.4739d-7,4.4869d-7,4.4933d-7,
     .4.4013d-7,4.5602d-7,4.4005d-7,5.1917d-7,5.0109d-7,4.9425d-7,
     .5.4399d-7,5.5298d-7,5.4768d-7,5.4768d-7,4.4240d-7,4.5286d-7,
     .4.4852d-7,4.1314d-7,4.1714d-7,4.2403d-7,3.1823d-7,3.1518d-7,
     .3.1732d-7,3.2355d-7,4.6249d-7,1.0035d-6,1.0527d-6,7.1976d-7,
     .1.0115d-6,1.0157d-6,6.8607d-7,6.9940d-7,7.0617d-7,5.7594d-7,
     .5.2900d-7,3.7412d-7,7.3528d-7,7.2422d-7,7.5358d-7,8.0152d-7,
     .7.0388d-7,5.6625d-7,5.1448d-7,6.7259d-7,6.7163d-7,5.1019d-7,
     .3.5045d-7,6.0376d-7,6.2746d-7,6.1980d-7,5.3402d-7,5.3402d-7,
     .6.9800d-7,7.0186d-7,1.0495d-6,1.0995d-6,7.4822d-7,8.6732d-7,
     .1.3526d-6,1.3487d-6,9.0438d-7,8.5835d-7,1.2710d-6,1.1108d-6,
     .7.2892d-7,1.0333d-6,1.0109d-6,5.2518d-7,4.7768d-7,7.2774d-7,
     .7.8360d-7,7.9655d-7,9.5685d-7,6.2447d-7,6.1596d-7,4.9642d-7,
     .1.0375d-6,1.0375d-6,9.2143d-7,8.2830d-7,5.0472d-7,5.1034d-7,
     .7.6203d-7,7.5472d-7,7.0006d-7,5.2249d-7,1.3037d-6,1.3019d-6,
     .7.2998d-7,5.3148d-7,5.1160d-7,9.7807d-7,9.5022d-7,5.2220d-7,
     .8.8391d-7,8.7904d-7,6.8228d-7,6.7948d-7,8.3622d-7,4.8253d-7,
     .8.3503d-7,8.3503d-7,4.9941d-7,6.7176d-7,6.6819d-7,4.6647d-7,
     .7.8900d-7,7.8689d-7,5.7918d-7,4.5911d-7,4.7977d-7,6.8763d-7,
     .8.1041d-7,8.3743d-7,8.8814d-7,5.1817d-7,5.5784d-7,8.0036d-7,
     .9.1155d-7,1.1788d-6,1.1609d-6,5.9684d-7,8.7403d-7,1.1661d-6,
     .1.1410d-6,6.0537d-7,6.0537d-7,8.7897d-7,1.2508d-6,1.2577d-6,
     .5.4430d-7,5.4430d-7,1.3325d-6,1.3325d-6,4.9817d-7,7.0704d-7,
     .9.5022d-7,9.6732d-7,9.7268d-7,7.1109d-7,4.7233d-7,6.6273d-7,
     .7.4651d-7,5.7020d-7,6.4517d-7,1.2251d-6,1.2319d-6,9.7092d-7,
     .1.0395d-6,1.3572d-6,1.3224d-6,8.4537d-7,9.6400d-7,1.0142d-6,
     .8.0711d-7,1.2415d-6,1.2761d-6,1.2623d-6,9.1678d-7,7.9384d-7,
     .7.7869d-7,8.8056d-7,9.7054d-7,2.2485d-6,2.1257d-6,1.5791d-6,
     .1.3095d-6,1.3130d-6,1.3777d-6,1.3760d-6,7.2153d-7,5.3687d-7,
     .8.9639d-7,1.4466d-6,1.4327d-6,7.1054d-7,6.7176d-7,1.2753d-6,
     .1.3492d-6,1.1561d-6,1.1264d-6,5.0549d-7,1.1120d-6,1.2835d-6,
     .1.3103d-6,1.3139d-6,1.2853d-6,1.0189d-6,9.6344d-7,1.3659d-6,
     .1.3548d-6,1.1292d-6,1.3320d-6,1.5912d-6,2.2098d-6,2.2098d-6,
     .8.9060d-7,1.2575d-6,1.2575d-6,6.8724d-7,7.0623d-7,7.1996d-7,
     .1.5292d-6,1.6198d-6,1.5847d-6,1.5523d-6,1.0453d-6/

      UPSILONMU=0d0

      IF(MX.LT.M(1).OR.MX.GT.M(N))RETURN

      DO I=2,N
       IF(MX.LT.M(I))THEN
        UPSILONMU=X(I-1)+(MX-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
        RETURN
       ENDIF
      ENDDO

      END


      DOUBLE PRECISION FUNCTION DEV(N,NN)
      IMPLICIT NONE
      DOUBLE PRECISION N,NN

      IF(N.EQ.NN)THEN
       DEV=0d0
      ELSE
       DEV=(NN-N)/(DABS(N)+DABS(NN))
      ENDIF

      END


      DOUBLE PRECISION FUNCTION DDSIN(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,D

      D=DSIN(X)
      IF(DABS(D).LT.1d-15)THEN
       DDSIN=0d0
      ELSE
       DDSIN=D
      ENDIF

      END

      DOUBLE PRECISION FUNCTION DDCOS(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,D

      D=DCOS(X)
      IF(DABS(D).LT.1d-15)THEN
       DDCOS=0d0
      ELSE
       DDCOS=D
      ENDIF

      END

