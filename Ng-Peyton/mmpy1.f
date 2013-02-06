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
C*************     MMPY1  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       LOOP UNROLLING: LEVEL 1
C
C   INPUT PARAMETERS -
C       M               -   NUMBER OF ROWS IN X AND IN Y.
C       N               -   NUMBER OF COLUMNS IN X AND NUMBER OF ROWS
C                           IN A.
C       Q               -   NUMBER OF COLUMNS IN A AND Y.
C       XPNT(*)         -   XPNT(J+1) POINTS ONE LOCATION BEYOND THE
C                           END OF THE J-TH COLUMN OF X.  XPNT IS ALSO
C                           USED TO ACCESS THE ROWS OF A.
C       X(*)            -   CONTAINS THE COLUMNS OF X AND THE ROWS OF A.
C       LDY             -   LENGTH OF FIRST COLUMN OF Y.
C
C   UPDATED PARAMETERS -
C       Y(*)            -   ON OUTPUT, Y = Y + AX.
C
C***********************************************************************
C
      SUBROUTINE  MMPY1  (  M     , N     , Q     , XPNT  , X     ,
     &                      Y     , LDY                             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER             LDY   , M     , N     , Q
        INTEGER             XPNT(*)
        DOUBLE PRECISION    X(*)          , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             I1
        INTEGER             IY    , IYLAST, IYSTRT, IYSTOP, LENY  ,
     &                      MM    , XCOL  , YCOL
        DOUBLE PRECISION    A1
C
C***********************************************************************
C
        MM = M
        IYLAST = 0
        LENY = LDY
C       ------------------------------------
C       TO COMPUTE EACH COLUMN YCOL OF Y ...
C       ------------------------------------
        DO  300  YCOL = 1, Q
            IYSTRT = IYLAST + 1
            IYSTOP = IYSTRT + MM - 1
            IYLAST = IYLAST + LENY
C           --------------------------------------------------
C           ... PERFORM THE APPROPRATE MATRIX VECTOR MULTIPLY:
C               X * A(*,YCOL).
C           --------------------------------------------------
            DO  200  XCOL = 1, N
                I1 = XPNT(XCOL+1) - MM
                A1  = - X(I1)
                DO  100  IY = IYSTRT, IYSTOP
                    Y(IY) = Y(IY) + A1 * X(I1)
                    I1 = I1 + 1
  100           CONTINUE
  200       CONTINUE
            MM = MM - 1
            LENY = LENY - 1
  300   CONTINUE
C
        RETURN
        END
