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
C**************     MMPY  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C   INPUT PARAMETERS -
C       M               -   NUMBER OF ROWS IN X AND IN Y.
C       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS
C                           IN A.
C       Q               -   NUMBER OF COLUMNS IN A AND Y.
C       SPLIT(*)        -   BLOCK PARTITIONING OF X.
C       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
C                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
C                           USED TO ACCESS THE ROWS OF A.
C       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
C       LDY             -   LENGTH OF FIRST COLUMN OF Y.
C       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY,
C                           WITH LEVEL N LOOP UNROLLING.
C
C   UPDATED PARAMETERS -
C       Y(*)            -   ON OUTPUT, Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  MMPY   (  M     , N     , Q     , SPLIT , XPNT  ,
     &                      X     , Y     , LDY   , MMPYN           )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        EXTERNAL            MMPYN
        INTEGER             LDY   , M     , N     , Q
        INTEGER             SPLIT(*)      , XPNT(*)
        DOUBLE PRECISION    X(*)          , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             BLK   , FSTCOL, NN
C
C***********************************************************************
C
        BLK = 1
        FSTCOL = 1
  100   CONTINUE
        IF  ( FSTCOL .LE. N )  THEN
            NN = SPLIT(BLK)
            CALL  MMPYN ( M, NN, Q, XPNT(FSTCOL), X, Y, LDY )
            FSTCOL = FSTCOL + NN
            BLK = BLK + 1
            GO TO 100
        ENDIF
        RETURN
C
      END
