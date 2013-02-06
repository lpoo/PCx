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
C
C     ------------------------------------
C     CONSTRUCT RIGHT HAND SIDE VECTOR ...
C     ------------------------------------
C
      SUBROUTINE  GETRHS (  N, XADJF, ADJF, ANZF, SOL, RHS )
C
        INTEGER             N
        INTEGER             XADJF(*), ADJF(*)
        DOUBLE PRECISION    ANZF(*), SOL(*), RHS(*)
C
        INTEGER             II, J
        DOUBLE PRECISION    T
C
C       ---------------
C       INITIALIZATION.
C       ---------------
        DO  100  J = 1, N
            RHS(J) = 0.0
  100   CONTINUE
C
        DO  300  J = 1, N
C           -------------------
C           FOR EACH COLUMN ...
C           -------------------
            T = SOL(J)
            DO  200  II = XADJF(J), XADJF(J+1)-1
                RHS(ADJF(II)) = RHS(ADJF(II)) + T*ANZF(II)
  200       CONTINUE
  300   CONTINUE
C
        RETURN
      END
