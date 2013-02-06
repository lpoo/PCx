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
C     -------------------------------------------
C     INSERT DIAGONAL NONZEROS INTO STRUCTURE AND
C     CREATE NUMERICAL VALUES.
C     -------------------------------------------
C
      SUBROUTINE  CREATE (  N, XADJ, ADJ, XADJF, ADJF, ANZF )
C
        INTEGER             N
        INTEGER             XADJ(*), ADJ(*)
        INTEGER             XADJF(*), ADJF(*)
        DOUBLE PRECISION    ANZF(*)
C
        INTEGER             II, J, JSTOP, JSTRT, NEXT, NOFNZ
C
        NEXT = 1
        XADJF(1) = NEXT
C
        JSTRT = XADJ(1)
        DO  200  J = 1, N
            JSTOP = XADJ(J+1) - 1
            NOFNZ = JSTOP - JSTRT + 1
C
            ADJF(NEXT) = J
            ANZF(NEXT) = NOFNZ + 1
            NEXT = NEXT + 1
C
            DO  100  II = JSTRT, JSTOP
                ADJF(NEXT) = ADJ(II)
                ANZF(NEXT) = -1.0
                NEXT = NEXT + 1
  100       CONTINUE
            XADJF(J+1) = NEXT
            JSTRT = JSTOP + 1
  200   CONTINUE
C
        RETURN
      END
