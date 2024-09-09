      SUBROUTINE LPTSOLVE(IFAIL)

      IMPLICIT NONE 

      INTEGER IFAIL,JM,JL,I,N,T
      PARAMETER (N=2)
      DOUBLE PRECISION INV(N,N),JAC(N,N),X(N),F(N),EPS,DET,DEV,TMAX
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU

      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL

      EPS=1d-8
      TMAX=100
      T=0
      JL=JL+1
!      WRITE(0,*)"CALL LPTSOLVE"
!      WRITE(0,*)""
!      WRITE(0,*)"T =",T
!      WRITE(0,*)""

      X(1)=LPPM
      X(2)=LTTM
!      WRITE(0,*)"X =",X
!      WRITE(0,*)""

      CALL LPTVAR(N,X,F,EPS,IFAIL)
      IF(IFAIL.NE.0)RETURN
!      WRITE(0,*)"F =",F

      DEV=DSQRT(F(1)**2+F(2)**2)
!      WRITE(0,*)"DEV =",DEV
!      WRITE(0,*)""

 1    T=T+1
!      WRITE(0,*)"T =",T
!      WRITE(0,*)""
      CALL LPTJAC(N,JAC,X,EPS,IFAIL)
      IF(IFAIL.NE.0)RETURN

      DET=JAC(1,1)*JAC(2,2)-JAC(1,2)*JAC(2,1)
      INV(1,1)=JAC(2,2)/DET
      INV(1,2)=-JAC(1,2)/DET
      INV(2,1)=-JAC(2,1)/DET
      INV(2,2)=JAC(1,1)/DET

      DO I=1,N
       X(I)=X(I)-INV(I,1)*F(1)-INV(I,2)*F(2)
      ENDDO
!      WRITE(0,*)"X =",X
!      WRITE(0,*)""

      CALL LPTVAR(N,X,F,EPS,IFAIL)
      IF(IFAIL.NE.0)RETURN
!      WRITE(0,*)"F =",F

      DEV=DSQRT(F(1)**2+F(2)**2)
!      WRITE(0,*)"DEV =",DEV
!      WRITE(0,*)""
      IF(DEV.GT.EPS .AND. T.LE.TMAX)THEN
       GOTO 1
      ELSEIF(T.GT.TMAX)THEN
       IFAIL=1
       RETURN
      ENDIF

      LPPMES=X(1)
      LPPT=LPPT+LPPMES
      LPPM=LPPT/JL
      LTTMES=X(2)
      LTTT=LTTT+LTTMES
      LTTM=LTTT/JL

      END


      SUBROUTINE LPTJAC(N,JAC,X,EPS,IFAIL)

      IMPLICIT NONE 

      INTEGER N,I,J,IFAIL
      DOUBLE PRECISION JAC(N,N),X(N),F(N),EPS,T,H

      DO J=1,N
       T=X(J)
       IF(T.EQ.0d0)THEN
        DO I=1,N
         JAC(I,J)=0d0
        ENDDO
       ELSE
        H=EPS*DABS(T)
        X(J)=T+H
        CALL LPTVAR(N,X,F,EPS,IFAIL)
        IF(IFAIL.NE.0)RETURN
        X(J)=T
        DO I=1,N
         IF(F(I).EQ.0d0)THEN
          JAC(I,J)=0d0
         ELSE
          JAC(I,J)=F(I)/(2d0*H)
         ENDIF
        ENDDO
        X(J)=T-H
        CALL LPTVAR(N,X,F,EPS,IFAIL)
        IF(IFAIL.NE.0)RETURN
        X(J)=T
        DO I=1,N
         IF(F(I).EQ.0d0)THEN
          JAC(I,J)=0d0
         ELSE
          JAC(I,J)=JAC(I,J)-F(I)/(2d0*H)
         ENDIF
        ENDDO
       ENDIF
      ENDDO

      END


      SUBROUTINE LPTVAR(N,X,F,EPS,IFAIL)

      IMPLICIT NONE 

      INTEGER N,IFAIL,JM,JL,GMFLAG,NN
      PARAMETER (NN=15)

      DOUBLE PRECISION X(N),F(N),EPS,X1,X2,Y(NN)
      DOUBLE PRECISION PI,COEF
      DOUBLE PRECISION G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      DOUBLE PRECISION MSUSYEFF,MMESS,N5,MGUT,MSREF,D,DMIN
      DOUBLE PRECISION LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      DOUBLE PRECISION LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      DOUBLE PRECISION MSM,MST,LM,LT,KM,KT,HTM,HTT,LPPM,LPPT,LTTM,LTTT

      COMMON/GMSCEN/MSREF,D,DMIN,GMFLAG
      COMMON/MESCOUP/G1MES,G2MES,G3MES,LMES,KMES,HTMES,HBMES,HLMES
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTGUT,HBGUT,HLGUT
      COMMON/MESGUT/LPPMES,LTTMES,LUMES,LDMES,LTMES,LBMES,LLMES,XIU
      COMMON/GUTMES/LPPGUT,LTTGUT,LUGUT,LDGUT,LTGUT,LBGUT,LLGUT
      COMMON/MESCAL/MSUSYEFF,MMESS,N5
      COMMON/MGUT/MGUT
      COMMON/GMSAVE/MSM,MST,LM,LT,KM,KT,HTM,HTT,
     . LPPM,LPPT,LTTM,LTTT,JM,JL

      EXTERNAL DERIVSGMGUT,RKQSGMGUT

      PI=4d0*DATAN(1d0)
      COEF=1d0/(16d0*PI**2)

* Definition of the couplings squared Y(I) at MMESS

      Y(1)=G1MES
      Y(2)=G2MES
      Y(3)=G3MES
      Y(4)=DSQRT(LM)
      Y(5)=KM
      Y(6)=DSQRT(HTM)
      Y(7)=DSQRT(HBMES)
      Y(8)=DSQRT(HLMES)
      Y(9)=X(1)
      Y(10)=X(2)
      Y(11)=LUMES
      Y(12)=LDMES
      Y(13)=LTMES
      Y(14)=LBMES
      Y(15)=LLMES

      X1=0d0
      X2=(3d0/G1MES-5d0/G2MES)/28d0

!      WRITE(0,*)"CALL LPTVAR"
!      WRITE(0,*)""
!      WRITE(0,*)"MMESS =",MMESS
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L =",Y(4)
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",Y(6)
!      WRITE(0,*)"HB =",Y(7)
!      WRITE(0,*)"HL =",Y(8)
!      WRITE(0,*)"LPP =",Y(9)
!      WRITE(0,*)"LTT =",Y(10)
!      WRITE(0,*)"LU =",Y(11)
!      WRITE(0,*)"LD =",Y(12)
!      WRITE(0,*)"LT =",Y(13)
!      WRITE(0,*)"LB =",Y(14)
!      WRITE(0,*)"LL =",Y(15)
!      WRITE(0,*)""

      CALL ODEINTGMGUT(Y,NN,X1,X2,EPS,DERIVSGMGUT,RKQSGMGUT,IFAIL)

      IF(IFAIL.GT.0)THEN
!       WRITE(0,*)"IFAIL =",IFAIL
!       WRITE(0,*)""
       IFAIL=2
       RETURN
      ELSE
       IFAIL=0
      ENDIF

* The GUT scale in GeV:

      MGUT=MMESS*DEXP(8d0*PI**2*X2)

* Couplings at the GUT scale

      G1GUT=Y(1)
      G2GUT=Y(2)
      G3GUT=Y(3)
      LGUT=Y(4)
      KGUT=Y(5)
      HTGUT=Y(6)
      HBGUT=Y(7)
      HLGUT=Y(8)
      LPPGUT=Y(9)
      LTTGUT=Y(10)
      LUGUT=Y(11)
      LDGUT=Y(12)
      LTGUT=Y(13)
      LBGUT=Y(14)
      LLGUT=Y(15)

      F(1)=LPPGUT-XIU
      F(2)=LTTGUT-XIU

!      WRITE(0,*)"MGUT =",MGUT
!      WRITE(0,*)"G1 =",5d0/3d0*Y(1)
!      WRITE(0,*)"G2 =",Y(2)
!      WRITE(0,*)"G3 =",Y(3)
!      WRITE(0,*)"L =",Y(4)
!      WRITE(0,*)"K =",Y(5)
!      WRITE(0,*)"HT =",Y(6)
!      WRITE(0,*)"HB =",Y(7)
!      WRITE(0,*)"HL =",Y(8)
!      WRITE(0,*)"LPP =",Y(9)
!      WRITE(0,*)"LTT =",Y(10)
!      WRITE(0,*)"LU =",Y(11)
!      WRITE(0,*)"LD =",Y(12)
!      WRITE(0,*)"LT =",Y(13)
!      WRITE(0,*)"LB =",Y(14)
!      WRITE(0,*)"LL =",Y(15)
!      WRITE(0,*)""

      END
