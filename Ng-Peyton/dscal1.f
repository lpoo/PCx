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
C******     DSCAL1 .... SCALE A VECTOR                     **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE COMPUTES A <-- AX, WHERE A IS A
C               SCALAR AND X IS A VECTOR.
C
C     INPUT PARAMETERS -
C        N - LENGTH OF THE VECTOR X.
C        A - SCALAR MULIPLIER.
C        X - VECTOR TO BE SCALED.
C
C     OUTPUT PARAMETERS - 
C        X - REPLACED BY THE SCALED VECTOR, AX.
C
C***********************************************************************
C
      SUBROUTINE  DSCAL1 ( N, A, X )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
      INTEGER             N
      DOUBLE PRECISION    A, X(N)
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
      INTEGER             I
C
C***********************************************************************
C
      DO  100  I = 1, N
          X(I) = A * X(I)
  100 CONTINUE
      RETURN
      END
