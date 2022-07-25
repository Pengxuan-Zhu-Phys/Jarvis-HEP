C Various CERNLIB mathematical routines needed for the optimum interval,
C as modified from 2006_src.tar.gz to not call other cernlib routines. 
C See http://cernlib.web.cern.ch/cernlib/download/2006_source/tar/
C
      FUNCTION DFREQ(X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
C
      DIMENSION P1(0:3),Q1(0:3),P2(0:7),Q2(0:7),P3(0:4),Q3(0:4)

      PARAMETER(Z1 = 1, HF = Z1/2)
      PARAMETER(C1 = 0.56418 95835 47756 29D0)
      PARAMETER(W2 = 1.41421 35623 73095 05D0, RW2 = 1/W2)

      DATA (P1(I),Q1(I),I=0,3)
     +/+2.42667 95523 05317 5D+2, +2.15058 87586 98612 0D+2,
     1 +2.19792 61618 29415 2D+1, +9.11649 05404 51490 1D+1,
     2 +6.99638 34886 19135 5D+0, +1.50827 97630 40778 7D+1,
     3 -3.56098 43701 81538 5D-2, +1/

      DATA (P2(I),Q2(I),I=0,7)
     +/+3.00459 26102 01616 01D+2, +3.00459 26095 69832 93D+2,
     1 +4.51918 95371 18729 42D+2, +7.90950 92532 78980 27D+2,
     2 +3.39320 81673 43436 87D+2, +9.31354 09485 06096 21D+2,
     3 +1.52989 28504 69404 04D+2, +6.38980 26446 56311 67D+2,
     4 +4.31622 27222 05673 53D+1, +2.77585 44474 39876 43D+2,
     5 +7.21175 82508 83093 66D+0, +7.70001 52935 22947 30D+1,
     6 +5.64195 51747 89739 71D-1, +1.27827 27319 62942 35D+1,
     7 -1.36864 85738 27167 07D-7, +1/

      DATA (P3(I),Q3(I),I=0,4)
     +/-2.99610 70770 35421 74D-3, +1.06209 23052 84679 18D-2,
     1 -4.94730 91062 32507 34D-2, +1.91308 92610 78298 41D-1,
     2 -2.26956 59353 96869 30D-1, +1.05167 51070 67932 07D+0,
     3 -2.78661 30860 96477 88D-1, +1.98733 20181 71352 56D+0,
     4 -2.23192 45973 41846 86D-2, +1/

      V=RW2*ABS(X)
      IF(V .LT. HF) THEN
       Y=V**2
       AP=P1(3)
       AQ=Q1(3)
       DO 1 I = 2,0,-1
       AP=P1(I)+Y*AP
    1  AQ=Q1(I)+Y*AQ
       H=V*AP/AQ
       HC=1-H
      ELSEIF(V .LT. 4) THEN
       AP=P2(7)
       AQ=Q2(7)
       DO 2 I = 6,0,-1
       AP=P2(I)+V*AP
    2  AQ=Q2(I)+V*AQ
       HC=EXP(-V**2)*AP/AQ
       H=1-HC
      ELSE
       Y=1/V**2
       AP=P3(4)
       AQ=Q3(4)
       DO 3 I = 3,0,-1
       AP=P3(I)+Y*AP
    3  AQ=Q3(I)+Y*AQ
       HC=EXP(-V**2)*(C1+Y*AP/AQ)/V
       H=1-HC
      ENDIF
      IF(X .GT. 0) THEN
       DFREQ=HF+HF*H
      ELSE
       DFREQ=HF*HC
      ENDIF
      RETURN
      END
      FUNCTION DGAUSN(P)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     Computes a "Normal Deviate"
C     Based on G.W. Hill & A.W. Davis, Algorithm 442 Normal Deviate
C     Collected Algorithms from CACM

C
      PARAMETER (C = 2.50662 82746 31000 50D0)
      PARAMETER (Z1 = 1, HF = Z1/2, C1 = 3*Z1/4, C2 = 7*Z1/8, C3 = Z1/3)

      IF(P .LE. 0 .OR. P .GE. 1) THEN
       H=0
       WRITE(6,101) P
      ELSEIF(P .EQ. HF) THEN
       H=0
      ELSE
       X=P
       IF(P .GT. HF) X=1-P
       X=SQRT(-2*LOG(X))
       X=X-((7.47395*X+494.877)*X+1637.720)/
     1     (((X+117.9407)*X+908.401)*X+659.935)
       IF(P .LT. HF) X=-X
       S=X**2
       Z=C*(P-DFREQ(X))*EXP(HF*S)
       H=(((((C1*S+C2)*Z+X)*X+HF)*C3*Z+HF*X)*Z+1)*Z+X
      ENDIF
      DGAUSN=H
      RETURN
  101 FORMAT('ARGUMENT P =',1P,D15.8,' NOT IN RANGE')
      END

      FUNCTION GAMDIS(X,A)
      IMPLICIT REAL (A-H,O-Z)
C     Calculates the gamma distribution function
C
C        G(x,a) = (1/gamma(a)) * int(0,x)[exp(-t) * t**(a-1)] dt.
C
C     Based on
C        W. Gautschi, ALGORITHM 542 Incomplete Gamma Functions,
C        ACM Trans. Math. Software 5 (1979) 482-489.

      DIMENSION C(14)

      PARAMETER (EPS = 1E-5, EPS1 = 5E-7)
      PARAMETER (ALH = -0.69314 72)
      PARAMETER (Z1 = 1, HALF = Z1/2, QUAR = Z1/4)
      PARAMETER (C1 = 3*Z1/2, KMAX = 300)

      DATA C
     1/ 0.5772157,-0.6558781,-0.0420026, 0.1665386,-0.0421977,
     2 -0.0096220, 0.0072189,-0.0011652,-0.0002152, 0.0001281,
     3 -0.0000201,-0.0000013, 0.0000011,-0.0000002/

      HST=0
      IF(X .EQ. 0) GO TO 99
      IF(X .LT. 0 .OR. A .LE. 0) THEN
       WRITE(6,101) X,A
       GO TO 99
      ELSE
       ALX=LOG(X)
      ENDIF
      IF(X .LT. QUAR) THEN
       ALFA=ALH/ALX
      ELSE
       ALFA=X+QUAR
      ENDIF
      IF(A .GT. ALFA) THEN
       TERM=1
       SUM=1
       DO 1 K = 1,KMAX
       TERM=X*TERM/(A+K)
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 2
    1  CONTINUE
       GO TO 98
    2  HST=SUM*EXP(A*ALX-X-ALOGAM(1+A))
      ELSEIF(X .GT. C1) THEN
       P=0
       S=1-A
       Q=(X+S)*(X-1-A)
       R=4*(X+S)
       TERM=1
       SUM=1
       RHO=0
       DO 3 K = 2,KMAX
       P=P+S
       Q=Q+R
       R=R+8
       S=S+2
       T=P*(1+RHO)
       RHO=T/(Q-T)
       TERM=RHO*TERM
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS*SUM) GO TO 4
    3  CONTINUE
       GO TO 98
    4  HST=1-(A*SUM/(X+1-A))*EXP(A*ALX-X-ALOGAM(1+A))
      ELSE
       IF(A .LT. HALF) THEN
        SUM=C(14)
        DO 12 K = 13,1,-1
   12   SUM=A*SUM+C(K)
        GA=-SUM/(1+A*SUM)
        Y=A*ALX
        IF(ABS(Y) .GE. 1) THEN
         U=GA-(EXP(Y)-1)/A
        ELSE
         SUM=1
         TERM=1
         DO 7 K = 2,KMAX
         TERM=Y*TERM/K
         SUM=SUM+TERM
         IF(ABS(TERM) .LE. EPS1*SUM) GO TO 8
    7    CONTINUE
         GO TO 98
    8    U=GA-SUM*ALX
        ENDIF
       ELSE
        U=GAMMA(A)-EXP(A*ALX)/A
       ENDIF
       P=A*X
       Q=A+1
       R=A+3
       TERM=1
       SUM=1
       DO 9 K = 2,KMAX
       P=P+X
       Q=Q+R
       R=R+2
       TERM=-P*TERM/Q
       SUM=SUM+TERM
       IF(ABS(TERM) .LE. EPS1*SUM) GO TO 10
    9  CONTINUE
       GO TO 98
   10  HST=1-A*(U+SUM*EXP((1+A)*ALX)/(1+A))/GAMMA(1+A)
      ENDIF
   99 GAMDIS=HST
      RETURN

   98 WRITE(6,102) X,A
      GO TO 99
  101 FORMAT('ILLEGAL ARGUMENT(S) X = ',E15.8,'  A = ',E15.8)
  102 FORMAT('PROBLEMS WITH CONVERGENCE, X = ',E15.8,'  A = ',E15.8)
      END
      SUBROUTINE RZERO(A,B,X0,R,EPS,MXF,F)
      IMPLICIT REAL (A-H,O-Z)

      EXTERNAL F
 
      PARAMETER (ONE = 1, HALF = ONE/2)
 
      XA=MIN(A,B)
      XB=MAX(A,B)
      FA=F(XA,1)
      FB=F(XB,2)
      IF(FA*FB .GT. 0) GO TO 5
      MC=0
 
    1 X0=HALF*(XA+XB)
      R=X0-XA
      EE=EPS*(ABS(X0)+1)
      IF(R .LE. EE) GO TO 4
      F1=FA
      X1=XA
      F2=FB
      X2=XB
 
    2 FX=F(X0,2)
      MC=MC+1
      IF(MC .GT. MXF) GO TO 6
      IF(FX*FA .GT. 0) THEN
       XA=X0
       FA=FX
      ELSE
       XB=X0
       FB=FX
      END IF
 
    3 U1=F1-F2
      U2=X1-X2
      U3=F2-FX
      U4=X2-X0
      IF(U2 .EQ. 0 .OR. U4 .EQ. 0) GO TO 1
      F3=FX
      X3=X0
      U1=U1/U2
      U2=U3/U4
      CA=U1-U2
      CB=(X1+X2)*U2-(X2+X0)*U1
      CC=(X1-X0)*F1-X1*(CA*X1+CB)
      IF(CA .EQ. 0) THEN
       IF(CB .EQ. 0) GO TO 1
       X0=-CC/CB
      ELSE
       U3=CB/(2*CA)
       U4=U3*U3-CC/CA
       IF(U4 .LT. 0) GO TO 1
       X0=-U3+SIGN(SQRT(U4),X0+U3)
      END IF
      IF(X0 .LT. XA .OR. X0 .GT. XB) GO TO 1
 
      R=MIN(ABS(X0-X3),ABS(X0-X2))
      EE=EPS*(ABS(X0)+1)
      IF(R .GT. EE) THEN
       F1=F2
       X1=X2
       F2=F3
       X2=X3
       GO TO 2
      END IF
 
      FX=F(X0,2)
      IF(FX .EQ. 0) GO TO 4
      IF(FX*FA .LT. 0) THEN
       XX=X0-EE
       IF(XX .LE. XA) GO TO 4
       FF=F(XX,2)
       FB=FF
       XB=XX
      ELSE
       XX=X0+EE
       IF(XX .GE. XB) GO TO 4
       FF=F(XX,2)
       FA=FF
       XA=XX
      END IF
      IF(FX*FF .GT. 0) THEN
       MC=MC+2
       IF(MC .GT. MXF) GO TO 6
       F1=F3
       X1=X3
       F2=FX
       X2=X0
       X0=XX
       FX=FF
       GO TO 3
      END IF
 
    4 R=EE
      FF=F(X0,3)
      RETURN
    5 Continue
      WRITE(*,100)
      R=-2*(XB-XA)
      X0=0
      RETURN
    6 Continue
      WRITE(*,101)
      R=-HALF*ABS(XB-XA)
      X0=0
      RETURN
  100 FORMAT(1X,'***** Modified RZERO ... F(A) AND F(B)',
     1          ' HAVE THE SAME SIGN')
  101 FORMAT(1X,'***** Modified RZERO ... TOO MANY FUNCTION CALLS')
      END

      FUNCTION ALGAMA(X)
C
C
      DIMENSION P1(7),Q1(7),P2(7),Q2(7),P3(7),Q3(7),C(5)

      PARAMETER (Z1 = 1, HF = Z1/2, HF1 = 1+HF)
      DATA P1
     1/+3.84287 36567 45991D+0, +5.27068 93753 00983D+1,
     2 +5.55840 45723 51531D+1, -2.15135 13573 72570D+2,
     3 -2.45872 61722 29242D+2, -5.75008 93603 04123D+1,
     4 -2.33590 98949 51284D+0/
      DATA Q1
     1/+1.00000 00000 00000D+0, +3.37330 47907 07074D+1,
     2 +1.93877 84034 37713D+2, +3.08829 54973 42428D+2,
     3 +1.50068 39064 89095D+2, +2.01068 51344 33395D+1,
     4 +4.57174 20282 50299D-1/
      DATA P2
     1/+4.87402 01396 83863 6D+0, +2.48845 25168 57407 6D+2,
     2 +2.17973 66058 89591 5D+3, +3.79751 24011 52511 8D+3,
     3 -1.97780 70769 84164 6D+3, -3.69298 34005 59128 2D+3,
     4 -5.60177 73537 80387 7D+2/
      DATA Q2
     1/+1.00000 00000 00000 0D+0, +9.50999 17418 20893 8D+1,
     2 +1.56120 45277 92863 5D+3, +7.23400 87928 94807 1D+3,
     3 +1.04595 76594 05895 9D+4, +4.16994 15153 20023 1D+3,
     4 +2.76785 83623 80410 1D+2/
      DATA P3
     1/-6.88062 40094 59425D+3, -4.30699 69819 57098D+5,
     2 -4.75045 94653 43956D+6, -2.94234 45930 32234D+6,
     3 +3.63218 04931 54257D+7, -3.35677 82814 54576D+6,
     4 -2.48043 69488 28593D+7/
      DATA Q3
     1/+1.00000 00000 00000D+0, -1.42168 29839 65146D+3,
     2 -1.55528 90280 85353D+5, -3.41525 17108 01107D+6,
     3 -2.09696 23255 80444D+7, -3.45441 75093 34395D+7,
     4 -9.16055 82863 71317D+6/
      DATA C
     1/ 1.12249 21356 561D-1,  7.95916 92961 204D-2,
     1 -1.70877 94611 020D-3,  9.18938 53320 467D-1,
     2  1.34699 05627 879D+0/


      ENTRY ALOGAM(X)

      IF(X .LE. 0) THEN
       H=0
       WRITE(6,101) X
      ELSE IF(X .EQ. 1 .OR. X .EQ. 2) THEN
       H=0
      ELSE IF(X .LE. HF) THEN
       Y=X+1
       AP=P1(1)
       AQ=Q1(1)
       DO 2 I = 2,7
       AP=P1(I)+Y*AP
    2  AQ=Q1(I)+Y*AQ
       H=-LOG(X)+X*AP/AQ
      ELSE IF(X .LE. HF1) THEN
       AP=P1(1)
       AQ=Q1(1)
       DO 3 I = 2,7
       AP=P1(I)+X*AP
    3  AQ=Q1(I)+X*AQ
       H=(X-1)*AP/AQ
      ELSE IF(X .LE. 4) THEN
       AP=P2(1)
       AQ=Q2(1)
       DO 4 I = 2,7
       AP=P2(I)+X*AP
    4  AQ=Q2(I)+X*AQ
       H=(X-2)*AP/AQ
      ELSE IF(X .LE. 12) THEN
       AP=P3(1)
       AQ=Q3(1)
       DO 5 I = 2,7
       AP=P3(I)+X*AP
    5  AQ=Q3(I)+X*AQ
       H=AP/AQ
      ELSE
       Y=1/X**2
       H=(X-HF)*LOG(X)-X+C(4)+(C(1)+Y*(C(2)+Y*C(3)))/
     1                                        ((C(5)+Y)*X)
      ENDIF
      ALGAMA=H
      RETURN
  101 FORMAT('NON-POSITIVE ARGUMENT  X = ',1P,E15.6)
      END
      FUNCTION GAMMA(X)
      IMPLICIT REAL (A-H,O-Z)
C
      DIMENSION C(0:15)

      DATA C( 0) /3.65738 77250 83382 44D0/
      DATA C( 1) /1.95754 34566 61268 27D0/
      DATA C( 2) /0.33829 71138 26160 39D0/
      DATA C( 3) /0.04208 95127 65575 49D0/
      DATA C( 4) /0.00428 76504 82129 09D0/
      DATA C( 5) /0.00036 52121 69294 62D0/
      DATA C( 6) /0.00002 74006 42226 42D0/
      DATA C( 7) /0.00000 18124 02333 65D0/
      DATA C( 8) /0.00000 01096 57758 66D0/
      DATA C( 9) /0.00000 00059 87184 05D0/
      DATA C(10) /0.00000 00003 07690 81D0/
      DATA C(11) /0.00000 00000 14317 93D0/
      DATA C(12) /0.00000 00000 00651 09D0/
      DATA C(13) /0.00000 00000 00025 96D0/
      DATA C(14) /0.00000 00000 00001 11D0/
      DATA C(15) /0.00000 00000 00000 04D0/

      U=X
      IF(U .LE. 0) THEN
       WRITE(6,101) U
       H=0
       B0=0
       B2=0
       F=0
       GO TO 9
      ENDIF
    8 F=1
      IF(U .LT. 3) THEN
       DO 1 I = 1,INT(4-U)
       F=F/U
    1  U=U+1
      ELSE
       DO 2 I = 1,INT(U-3)
       U=U-1
    2  F=F*U
      END IF
      H=U+U-7
      ALFA=H+H
      B1=0
      B2=0
      DO 3 I = 15,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    3 B1=B0
    9 GAMMA =F*(B0-H*B2)
      RETURN
  101 FORMAT('ARGUMENT IS NEGATIVE = ',1P,E15.1)
      END
      SUBROUTINE RANLUX(RVEC,LENV)
C         Subtract-and-borrow random number generator proposed by
C         Marsaglia and Zaman, implemented by F. James with the name
C         RCARRY in 1991, and later improved by Martin Luescher
C         in 1993 to produce "Luxury Pseudorandom Numbers".
C     Fortran 77 coded by F. James, 1993
C
C   LUXURY LEVELS.
C   ------ ------      The available luxury levels are:
C
C  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
C           and Zaman, very long period, but fails many tests.
C  level 1  (p=48): considerable improvement in quality over level 0,
C           now passes the gap test, but still fails spectral test.
C  level 2  (p=97): passes all known tests, but theoretically still
C           defective.
C  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
C           correlations have very small chance of being observed.
C  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
C
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RANLUX:                                  ++
C!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
C!!!                   32-bit random floating point numbers between  ++
C!!!                   zero (not included) and one (also not incl.). ++
C!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
C!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
C!!!               which is integer between zero and MAXLEV, or if   ++
C!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
C!!!               should be set to zero unless restarting at a break++ 
C!!!               point given by output of RLUXAT (see RLUXAT).     ++
C!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
C!!!               which can be used to restart the RANLUX generator ++
C!!!               at the current point by calling RLUXGO.  K1 and K2++
C!!!               specify how many numbers were generated since the ++
C!!!               initialization with LUX and INT.  The restarting  ++
C!!!               skips over  K1+K2*E9   numbers, so it can be long.++
C!!!   A more efficient but less convenient way of restarting is by: ++
C!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
C!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
C!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
C!!!                 32-bit integer seeds, to be used for restarting ++
C!!!      ISVEC must be dimensioned 25 in the calling program        ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
C                               default
C  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
Corresponds to p=24    48    97   223   389
C     time factor 1     2     3     6    10   on slow workstation
C                 1    1.5    2     3     5   on fast mainframe
C
C  NOTYET is .TRUE. if no initialization has been performed yet.
C              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
            TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
C
C          The Generator proper: "Subtract-with-borrow",
C          as proposed by Marsaglia and Zaman,
C          Florida State University, March, 1989
C
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
C  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
C        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
C        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
C
C           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
         NOTYET = .FALSE.
         TWOM24 = 1.
         DO 195 I= 1, 24
         NEXT(I) = I-1
  195    TWOM24 = TWOM24 * 0.5
         NEXT(1) = 24
         TWOM12 = TWOM24 * 4096.
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
C
C                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
C
C                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
C
C                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
      ELSE
          NSKIP = LUXLEV - 24
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +      JSEED, K1,K2
      ELSE
        JSEED = JSDFLT
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
C        If restarting at a break point, skip K1 + IGIGA*K2
C        Note that this is the number of numbers delivered to
C        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
C         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
C       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
