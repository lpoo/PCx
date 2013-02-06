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
C*************     MMPYI  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       MATRIX X HAS ONLY 1 COLUMN.
C
C   INPUT PARAMETERS -
C       M               -   NUMBER OF ROWS IN X AND IN Y.
C       Q               -   NUMBER OF COLUMNS IN A AND Y.
C       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
C                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
C                           USED TO ACCESS THE ROWS OF A.
C       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
C       IY(*)           -   IY(COL) POINTS TO THE BEGINNING OF COLUMN
C       RELIND(*)       -   RELATIVE INDICES.
C
C   UPDATED PARAMETERS -
C       Y(*)            -   ON OUTPUT, Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  MMPYI  (  M     , Q     , XPNT  , X     , IY    ,
     &                      Y     , RELIND                          )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER             M     , Q
        INTEGER             IY(*)         , RELIND(*)     ,
     &                      XPNT(*)
        DOUBLE PRECISION    X(*)      , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             COL   , I     , ISUB  , K     , YLAST
        DOUBLE PRECISION    A
C
C***********************************************************************
C
        DO  200  K = 1, Q
            COL = XPNT(K)
            YLAST = IY(COL+1) - 1
            A = - X(K)
CDIR$   IVDEP
            DO  100  I = K, M
                ISUB = XPNT(I)
                ISUB = YLAST - RELIND(ISUB)
                Y(ISUB) = Y(ISUB) + A*X(I)
  100       CONTINUE
  200   CONTINUE
        RETURN
C
      END
