C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C******     SMXPY1 .... MATRIX-VECTOR MULTIPLY            **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS A MATRIX-VECTOR MULTIPLY,
C               Y = Y + AX, ASSUMING DATA STRUCTURES USED IN
C               RECENTLY DEVELOPED SPARSE CHOLESKY CODES.  THE 
C               '1' SIGNIFIES NO LOOP UNROLLING, I.E., 
C               LOOP-UNROLLING TO LEVEL 1.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS.
C        N      - NUMBER OF COLUMNS.
C        Y      - M-VECTOR TO WHICH AX WILL BE ADDED.
C        APNT   - INDEX VECTOR FOR A.  XA(I) POINTS TO THE
C                 FIRST NONZERO IN COLUMN I OF A.
C        Y      - ON OUTPUT, CONTAINS Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  SMXPY1 ( M, N, Y, APNT, A )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
C
      INTEGER             M, N
C
      INTEGER             APNT(N)
C
      DOUBLE PRECISION    Y(M), A(*)
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
C
      INTEGER             I, II, J
C
      DOUBLE PRECISION    AMULT
C
C***********************************************************************
C
      DO  200  J = 1, N
          II = APNT(J+1) - M
          AMULT  = - A(II)
          DO  100  I = 1, M
              Y(I) = Y(I) + AMULT * A(II)
              II = II + 1
  100     CONTINUE
  200 CONTINUE
      RETURN
      END
