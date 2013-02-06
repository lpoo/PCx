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
C****     ORDNAT ..... NATURAL ORDERING                     ************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE RECORDS THE INITIAL ORDERING IN THE
C               ORDERING VECTORS PERM AND INVP.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE "NATURAL" ORDERING; I.E., THE INITIAL
C                 ORDERING.
C        INVP   - THE INVERSE OF PERM.
C
C***********************************************************************
C
      SUBROUTINE ORDNAT  (  NEQNS , PERM  , INVP    )
C
C***********************************************************************
C
        INTEGER    INVP(1) ,    PERM(1)
        INTEGER    NEQNS
C
        INTEGER     I
C
C***********************************************************************
C
        DO  700  I = 1, NEQNS
            PERM(I) = I
            INVP(I) = I
  700   CONTINUE
        RETURN
C
      END
