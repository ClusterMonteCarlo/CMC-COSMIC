C     ALGORITHM 642 COLLECTED ALGORITHMS FROM ACM.
C     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.12, NO. 2,
C     JUN., 1986, P. 150.
C   SUBROUTINE NAME     - CUBGCV
C
C--------------------------------------------------------------------------
C
C   COMPUTER            - VAX/DOUBLE
C
C   AUTHOR              - M.F.HUTCHINSON
C                         CSIRO DIVISION OF MATHEMATICS AND STATISTICS
C                         P.O. BOX 1965
C                         CANBERRA, ACT 2601
C                         AUSTRALIA
C
C   LATEST REVISION     - 15 AUGUST 1985
C
C   PURPOSE             - CUBIC SPLINE DATA SMOOTHER
C
C   USAGE               - CALL CUBGCV (X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C
C   ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE
C                           ABSCISSAE OF THE N DATA POINTS
C                           (X(I),F(I)) I=1..N. (INPUT) X
C                           MUST BE ORDERED SO THAT
C                           X(I) .LT. X(I+1).
C                F      - VECTOR OF LENGTH N CONTAINING THE
C                           ORDINATES (OR FUNCTION VALUES)
C                           OF THE N DATA POINTS (INPUT).
C                DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT)
C                           DF(I) IS THE RELATIVE STANDARD DEVIATION
C                           OF THE ERROR ASSOCIATED WITH DATA POINT I.
C                           EACH DF(I) MUST BE POSITIVE.  THE VALUES IN
C                           DF ARE SCALED BY THE SUBROUTINE SO THAT
C                           THEIR MEAN SQUARE VALUE IS 1, AND UNSCALED
C                           AGAIN ON NORMAL EXIT.
C                           THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED
C                           IN WK(7) ON NORMAL EXIT.
C                           IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN,
C                           THESE SHOULD BE PROVIDED IN DF AND THE ERROR
C                           VARIANCE PARAMETER VAR (SEE BELOW) SHOULD THEN
C                           BE SET TO 1.
C                           IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN,
C                           SET EACH DF(I)=1.
C                N      - NUMBER OF DATA POINTS (INPUT).
C                           N MUST BE .GE. 3.
C                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y
C                           IS A VECTOR OF LENGTH N. C IS
C                           AN N-1 BY 3 MATRIX. THE VALUE
C                           OF THE SPLINE APPROXIMATION AT T IS
C                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I)
C                           WHERE X(I).LE.T.LT.X(I+1) AND
C                           D = T-X(I).
C                IC     - ROW DIMENSION OF MATRIX C EXACTLY
C                           AS SPECIFIED IN THE DIMENSION
C                           STATEMENT IN THE CALLING PROGRAM. (INPUT)
C                VAR    - ERROR VARIANCE. (INPUT/OUTPUT)
C                           IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN
C                           THE SMOOTHING PARAMETER IS DETERMINED
C                           BY MINIMIZING THE GENERALIZED CROSS VALIDATION
C                           AND AN ESTIMATE OF THE ERROR VARIANCE IS
C                           RETURNED IN VAR.
C                           IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE
C                           SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE
C                           AN ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE
C                           MEAN SQUARE ERROR, AND VAR IS UNCHANGED.
C                           IN PARTICULAR, IF VAR IS ZERO, THEN AN
C                           INTERPOLATING NATURAL CUBIC SPLINE IS CALCULATED.
C                           VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD
C                           DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE).
C                JOB    - JOB SELECTION PARAMETER. (INPUT)
C                         JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR
C                           ESTIMATES ARE NOT REQUIRED IN SE.
C                         JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR
C                           ESTIMATES ARE REQUIRED IN SE.
C                SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD
C                           ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y.
C                           SE IS NOT REFERENCED IF JOB=0. (OUTPUT)
C                WK     - WORK VECTOR OF LENGTH 7*(N + 2). ON NORMAL EXIT THE
C                           FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:-
C
C                           WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1))
C                           WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF
C                                   FREEDOM OF THE RESIDUAL SUM OF SQUARES
C                           WK(3) = GENERALIZED CROSS VALIDATION
C                           WK(4) = MEAN SQUARE RESIDUAL
C                           WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR
C                                   AT THE DATA POINTS
C                           WK(6) = ESTIMATE OF THE ERROR VARIANCE
C                           WK(7) = MEAN SQUARE VALUE OF THE DF(I)
C
C                           IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC
C                           SPLINE HAS BEEN CALCULATED.
C                           IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES
C                           REGRESSION LINE HAS BEEN CALCULATED.
C                           WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF
C                           FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE
C                           USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION
C                           LINE IS CALCULATED.
C                           WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I)
C                           SCALED TO HAVE MEAN SQUARE VALUE 1.  THE
C                           UNSCALED VALUES OF WK(3),WK(4),WK(5) MAY BE
C                           CALCULATED BY DIVIDING BY WK(7).
C                           WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF
C                           VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH
C                           THE UNSCALED VALUES OF THE DF(I) TO FACILITATE
C                           COMPARISONS WITH A PRIORI VARIANCE ESTIMATES.
C
C                IER    - ERROR PARAMETER. (OUTPUT)
C                         TERMINAL ERROR
C                           IER = 129, IC IS LESS THAN N-1.
C                           IER = 130, N IS LESS THAN 3.
C                           IER = 131, INPUT ABSCISSAE ARE NOT
C                             ORDERED SO THAT X(I).LT.X(I+1).
C                           IER = 132, DF(I) IS NOT POSITIVE FOR SOME I.
C                           IER = 133, JOB IS NOT 0 OR 1.
C
C   PRECISION/HARDWARE  - DOUBLE
C
C   REQUIRED ROUTINES   - SPINT1,SPFIT1,SPCOF1,SPERR1
C
C   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE
C                SUBROUTINE IS PROPORTIONAL TO N.  THE SUBROUTINE
C                USES AN ALGORITHM DEVELOPED BY M.F. HUTCHINSON AND
C                F.R. DE HOOG, 'SMOOTHING NOISY DATA WITH SPLINE
C                FUNCTIONS', NUMER. MATH. (IN PRESS)
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CUBGCV(X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER)
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,JOB,IER
      DOUBLE PRECISION X(N),F(N),DF(N),Y(N),C(IC,3),SE(N),VAR,
     .                 WK(0:N+1,7)
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      DOUBLE PRECISION DELTA,ERR,GF1,GF2,GF3,GF4,R1,R2,R3,R4,TAU,RATIO,
     .                 AVH,AVDF,AVAR,ZERO,ONE,STAT(6),P,Q
C
      DATA RATIO/2.0D0/
      DATA TAU/1.618033989D0/
      DATA ZERO,ONE/0.0D0,1.0D0/
C
C---INITIALIZE---
      IER = 133
      IF (JOB.LT.0 .OR. JOB.GT.1) GO TO 140
      CALL SPINT1(X,AVH,F,DF,AVDF,N,Y,C,IC,WK,WK(0,4),IER)
      IF (IER.NE.0) GO TO 140
      AVAR = VAR
      IF (VAR.GT.ZERO) AVAR = VAR*AVDF*AVDF
C
C---CHECK FOR ZERO VARIANCE---
      IF (VAR.NE.ZERO) GO TO 10
      R1 = ZERO
      GO TO 90
C
C---FIND LOCAL MINIMUM OF GCV OR THE EXPECTED MEAN SQUARE ERROR---
   10 R1 = ONE
      R2 = RATIO*R1
      CALL SPFIT1(X,AVH,DF,N,R2,P,Q,GF2,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
   20 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      IF (GF1.GT.GF2) GO TO 30
C
C---EXIT IF P ZERO---
      IF (P.LE.ZERO) GO TO 100
      R2 = R1
      GF2 = GF1
      R1 = R1/RATIO
      GO TO 20

   30 R3 = RATIO*R2
   40 CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      IF (GF3.GT.GF2) GO TO 50
C
C---EXIT IF Q ZERO---
      IF (Q.LE.ZERO) GO TO 100
      R2 = R3
      GF2 = GF3
      R3 = RATIO*R3
      GO TO 40

   50 R2 = R3
      GF2 = GF3
      DELTA = (R2-R1)/TAU
      R4 = R1 + DELTA
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
C
C---GOLDEN SECTION SEARCH FOR LOCAL MINIMUM---
   60 IF (GF3.GT.GF4) GO TO 70
      R2 = R4
      GF2 = GF4
      R4 = R3
      GF4 = GF3
      DELTA = DELTA/TAU
      R3 = R2 - DELTA
      CALL SPFIT1(X,AVH,DF,N,R3,P,Q,GF3,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
      GO TO 80

   70 R1 = R3
      GF1 = GF3
      R3 = R4
      GF3 = GF4
      DELTA = DELTA/TAU
      R4 = R1 + DELTA
      CALL SPFIT1(X,AVH,DF,N,R4,P,Q,GF4,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
   80 ERR = (R2-R1)/ (R1+R2)
      IF (ERR*ERR+ONE.GT.ONE .AND. ERR.GT.1.0D-6) GO TO 60
      R1 = (R1+R2)*0.5D0
C
C---CALCULATE SPLINE COEFFICIENTS---
   90 CALL SPFIT1(X,AVH,DF,N,R1,P,Q,GF1,AVAR,STAT,Y,C,IC,WK,WK(0,4),
     .            WK(0,6),WK(0,7))
  100 CALL SPCOF1(X,AVH,F,DF,N,P,Q,Y,C,IC,WK(0,6),WK(0,7))
C
C---OPTIONALLY CALCULATE STANDARD ERROR ESTIMATES---
      IF (VAR.GE.ZERO) GO TO 110
      AVAR = STAT(6)
      VAR = AVAR/ (AVDF*AVDF)
  110 IF (JOB.EQ.1) CALL SPERR1(X,AVH,DF,N,WK,P,AVAR,SE)
C
C---UNSCALE DF---
      DO 120 I = 1,N
         DF(I) = DF(I)*AVDF
  120 CONTINUE
C
C--PUT STATISTICS IN WK---
      DO 130 I = 0,5
         WK(I,1) = STAT(I+1)
  130 CONTINUE
      WK(5,1) = STAT(6)/ (AVDF*AVDF)
      WK(6,1) = AVDF*AVDF
      GO TO 150
C
C---CHECK FOR ERROR CONDITION---
  140 CONTINUE
C     IF (IER.NE.0) CONTINUE
  150 RETURN
      END
      SUBROUTINE SPINT1(X,AVH,Y,DY,AVDY,N,A,C,IC,R,T,IER)
C
C INITIALIZES THE ARRAYS C, R AND T FOR ONE DIMENSIONAL CUBIC
C SMOOTHING SPLINE FITTING BY SUBROUTINE SPFIT1.  THE VALUES
C DF(I) ARE SCALED SO THAT THE SUM OF THEIR SQUARES IS N
C AND THE AVERAGE OF THE DIFFERENCES X(I+1) - X(I) IS CALCULATED
C IN AVH IN ORDER TO AVOID UNDERFLOW AND OVERFLOW PROBLEMS IN
C SPFIT1.
C
C SUBROUTINE SETS IER IF ELEMENTS OF X ARE NON-INCREASING,
C IF N IS LESS THAN 3, IF IC IS LESS THAN N-1 OR IF DY(I) IS
C NOT POSITIVE FOR SOME I.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N,IC,IER
      DOUBLE PRECISION X(N),Y(N),DY(N),A(N),C(IC,3),R(0:N+1,3),
     .                 T(0:N+1,2),AVH,AVDY
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION E,F,G,H,ZERO
      DATA ZERO/0.0D0/
C
C---INITIALIZATION AND INPUT CHECKING---
      IER = 0
      IF (N.LT.3) GO TO 60
      IF (IC.LT.N-1) GO TO 70
C
C---GET AVERAGE X SPACING IN AVH---
      G = ZERO
      DO 10 I = 1,N - 1
         H = X(I+1) - X(I)
         IF (H.LE.ZERO) GO TO 80
         G = G + H
   10 CONTINUE
      AVH = G/ (N-1)
C
C---SCALE RELATIVE WEIGHTS---
      G = ZERO
      DO 20 I = 1,N
         IF (DY(I).LE.ZERO) GO TO 90
         G = G + DY(I)*DY(I)
   20 CONTINUE
      AVDY = DSQRT(G/N)
C
      DO 30 I = 1,N
         DY(I) = DY(I)/AVDY
   30 CONTINUE
C
C---INITIALIZE H,F---
      H = (X(2)-X(1))/AVH
      F = (Y(2)-Y(1))/H
C
C---CALCULATE A,T,R---
      DO 40 I = 2,N - 1
         G = H
         H = (X(I+1)-X(I))/AVH
         E = F
         F = (Y(I+1)-Y(I))/H
         A(I) = F - E
         T(I,1) = 2.0D0* (G+H)/3.0D0
         T(I,2) = H/3.0D0
         R(I,3) = DY(I-1)/G
         R(I,1) = DY(I+1)/H
         R(I,2) = -DY(I)/G - DY(I)/H
   40 CONTINUE
C
C---CALCULATE C = R'*R---
      R(N,2) = ZERO
      R(N,3) = ZERO
      R(N+1,3) = ZERO
      DO 50 I = 2,N - 1
         C(I,1) = R(I,1)*R(I,1) + R(I,2)*R(I,2) + R(I,3)*R(I,3)
         C(I,2) = R(I,1)*R(I+1,2) + R(I,2)*R(I+1,3)
         C(I,3) = R(I,1)*R(I+2,3)
   50 CONTINUE
      RETURN
C
C---ERROR CONDITIONS---
   60 IER = 130
      RETURN

   70 IER = 129
      RETURN

   80 IER = 131
      RETURN

   90 IER = 132
      RETURN
      END
      SUBROUTINE SPFIT1(X,AVH,DY,N,RHO,P,Q,FUN,VAR,STAT,A,C,IC,R,T,U,V)
C
C FITS A CUBIC SMOOTHING SPLINE TO DATA WITH RELATIVE
C WEIGHTING DY FOR A GIVEN VALUE OF THE SMOOTHING PARAMETER
C RHO USING AN ALGORITHM BASED ON THAT OF C.H. REINSCH (1967),
C NUMER. MATH. 10, 177-183.
C
C THE TRACE OF THE INFLUENCE MATRIX IS CALCULATED USING AN
C ALGORITHM DEVELOPED BY M.F.HUTCHINSON AND F.R.DE HOOG (NUMER.
C MATH., IN PRESS), ENABLING THE GENERALIZED CROSS VALIDATION
C AND RELATED STATISTICS TO BE CALCULATED IN ORDER N OPERATIONS.
C
C THE ARRAYS A, C, R AND T ARE ASSUMED TO HAVE BEEN INITIALIZED
C BY THE SUBROUTINE SPINT1.  OVERFLOW AND UNDERFLOW PROBLEMS ARE
C AVOIDED BY USING P=RHO/(1 + RHO) AND Q=1/(1 + RHO) INSTEAD OF
C RHO AND BY SCALING THE DIFFERENCES X(I+1) - X(I) BY AVH.
C
C THE VALUES IN DF ARE ASSUMED TO HAVE BEEN SCALED SO THAT THE
C SUM OF THEIR SQUARED VALUES IS N.  THE VALUE IN VAR, WHEN IT IS
C NON-NEGATIVE, IS ASSUMED TO HAVE BEEN SCALED TO COMPENSATE FOR
C THE SCALING OF THE VALUES IN DF.
C
C THE VALUE RETURNED IN FUN IS AN ESTIMATE OF THE TRUE MEAN SQUARE
C WHEN VAR IS NON-NEGATIVE, AND IS THE GENERALIZED CROSS VALIDATION
C WHEN VAR IS NEGATIVE.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      DOUBLE PRECISION X(N),DY(N),RHO,STAT(6),A(N),C(IC,3),R(0:N+1,3),
     .                 T(0:N+1,2),U(0:N+1),V(0:N+1),FUN,VAR,AVH,P,Q
C
C---LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION E,F,G,H,ZERO,ONE,TWO,RHO1
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C
C---USE P AND Q INSTEAD OF RHO TO PREVENT OVERFLOW OR UNDERFLOW---
      RHO1 = ONE + RHO
      P = RHO/RHO1
      Q = ONE/RHO1
      IF (RHO1.EQ.ONE) P = ZERO
      IF (RHO1.EQ.RHO) Q = ZERO
C
C---RATIONAL CHOLESKY DECOMPOSITION OF P*C + Q*T---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 10 I = 0,1
         R(I,1) = ZERO
   10 CONTINUE
      DO 20 I = 2,N - 1
         R(I-2,3) = G*R(I-2,1)
         R(I-1,2) = F*R(I-1,1)
         R(I,1) = ONE/ (P*C(I,1)+Q*T(I,1)-F*R(I-1,2)-G*R(I-2,3))
         F = P*C(I,2) + Q*T(I,2) - H*R(I-1,2)
         G = H
         H = P*C(I,3)
   20 CONTINUE
C
C---SOLVE FOR U---
      U(0) = ZERO
      U(1) = ZERO
      DO 30 I = 2,N - 1
         U(I) = A(I) - R(I-1,2)*U(I-1) - R(I-2,3)*U(I-2)
   30 CONTINUE
      U(N) = ZERO
      U(N+1) = ZERO
      DO 40 I = N - 1,2,-1
         U(I) = R(I,1)*U(I) - R(I,2)*U(I+1) - R(I,3)*U(I+2)
   40 CONTINUE
C
C---CALCULATE RESIDUAL VECTOR V---
      E = ZERO
      H = ZERO
      DO 50 I = 1,N - 1
         G = H
         H = (U(I+1)-U(I))/ ((X(I+1)-X(I))/AVH)
         V(I) = DY(I)* (H-G)
         E = E + V(I)*V(I)
   50 CONTINUE
      V(N) = DY(N)* (-H)
      E = E + V(N)*V(N)
C
C---CALCULATE UPPER THREE BANDS OF INVERSE MATRIX---
      R(N,1) = ZERO
      R(N,2) = ZERO
      R(N+1,1) = ZERO
      DO 60 I = N - 1,2,-1
         G = R(I,2)
         H = R(I,3)
         R(I,2) = -G*R(I+1,1) - H*R(I+1,2)
         R(I,3) = -G*R(I+1,2) - H*R(I+2,1)
         R(I,1) = R(I,1) - G*R(I,2) - H*R(I,3)
   60 CONTINUE
C
C---CALCULATE TRACE---
      F = ZERO
      G = ZERO
      H = ZERO
      DO 70 I = 2,N - 1
         F = F + R(I,1)*C(I,1)
         G = G + R(I,2)*C(I,2)
         H = H + R(I,3)*C(I,3)
   70 CONTINUE
      F = F + TWO* (G+H)
C
C---CALCULATE STATISTICS---
      STAT(1) = P
      STAT(2) = F*P
      STAT(3) = N*E/ (F*F)
      STAT(4) = E*P*P/N
      STAT(6) = E*P/F
      IF (VAR.GE.ZERO) GO TO 80
      STAT(5) = STAT(6) - STAT(4)
      FUN = STAT(3)
      GO TO 90

   80 STAT(5) = DMAX1(STAT(4)-TWO*VAR*STAT(2)/N+VAR,ZERO)
      FUN = STAT(5)
   90 RETURN
      END
      SUBROUTINE SPERR1(X,AVH,DY,N,R,P,VAR,SE)
C
C CALCULATES BAYESIAN ESTIMATES OF THE STANDARD ERRORS OF THE FITTED
C VALUES OF A CUBIC SMOOTHING SPLINE BY CALCULATING THE DIAGONAL ELEMENTS
C OF THE INFLUENCE MATRIX.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER N
      DOUBLE PRECISION X(N),DY(N),R(0:N+1,3),SE(N),AVH,P,VAR
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION F,G,H,F1,G1,H1,ZERO,ONE
      DATA ZERO,ONE/0.0D0,1.0D0/
C
C---INITIALIZE---
      H = AVH/ (X(2)-X(1))
      SE(1) = ONE - P*DY(1)*DY(1)*H*H*R(2,1)
      R(1,1) = ZERO
      R(1,2) = ZERO
      R(1,3) = ZERO
C
C---CALCULATE DIAGONAL ELEMENTS---
      DO 10 I = 2,N - 1
         F = H
         H = AVH/ (X(I+1)-X(I))
         G = -F - H
         F1 = F*R(I-1,1) + G*R(I-1,2) + H*R(I-1,3)
         G1 = F*R(I-1,2) + G*R(I,1) + H*R(I,2)
         H1 = F*R(I-1,3) + G*R(I,2) + H*R(I+1,1)
         SE(I) = ONE - P*DY(I)*DY(I)* (F*F1+G*G1+H*H1)
   10 CONTINUE
      SE(N) = ONE - P*DY(N)*DY(N)*H*H*R(N-1,1)
C
C---CALCULATE STANDARD ERROR ESTIMATES---
      DO 20 I = 1,N
         SE(I) = DSQRT(DMAX1(SE(I)*VAR,ZERO))*DY(I)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE SPCOF1(X,AVH,Y,DY,N,P,Q,A,C,IC,U,V)
C
C CALCULATES COEFFICIENTS OF A CUBIC SMOOTHING SPLINE FROM
C PARAMETERS CALCULATED BY SUBROUTINE SPFIT1.
C
C---SPECIFICATIONS FOR ARGUMENTS---
      INTEGER IC,N
      DOUBLE PRECISION X(N),Y(N),DY(N),P,Q,A(N),C(IC,3),U(0:N+1),
     .                 V(0:N+1),AVH
C
C---SPECIFICATIONS FOR LOCAL VARIABLES---
      INTEGER I
      DOUBLE PRECISION H,QH
C
C---CALCULATE A---
      QH = Q/ (AVH*AVH)
      DO 10 I = 1,N
         A(I) = Y(I) - P*DY(I)*V(I)
         U(I) = QH*U(I)
   10 CONTINUE
C
C---CALCULATE C---
      DO 20 I = 1,N - 1
         H = X(I+1) - X(I)
         C(I,3) = (U(I+1)-U(I))/ (3.0D0*H)
         C(I,1) = (A(I+1)-A(I))/H - (H*C(I,3)+U(I))*H
         C(I,2) = U(I)
   20 CONTINUE
      RETURN
      END
