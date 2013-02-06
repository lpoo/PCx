C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  May 26, 1995
C   Authors:        Esmond G. Ng, Barry W. Peyton, and Guodong Zhang
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C*************     MMPY2  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       LOOP UNROLLING: LEVEL 2 UPDATING TWO COLUMNS AT A TIME
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
        SUBROUTINE  MMPY2  (  M     , N     , Q     , XPNT  , X     ,
     &                        Y     , LDY                             )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER               LDY   , M     , N     , Q
        INTEGER               XPNT(*)
        DOUBLE PRECISION      X(*)          , Y(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER               I     , J     , K     , QQ
        INTEGER               I1    , I2
        INTEGER               IYBEG , IYBEG1, IYBEG2, LENY  , MM    
        DOUBLE PRECISION      A1    , A2    , A9    , A10
        DOUBLE PRECISION      B1    , B2    , Y1    , Y2
C
C***********************************************************************
C
C       ----------------------------------------------------
C       COMPUTE EACH DIAGONAL ENTRY OF THE ODD COLUMNS OF Y.
C       ----------------------------------------------------
C
        MM = M
        QQ = MIN(M,Q)
        IYBEG = 1
        LENY = LDY - 1
        DO  200 J = 1, QQ-1 , 2
CDIR$   IVDEP
            DO  100  I = 1, N
                I1 = XPNT(I+1) - MM
                A1 = X(I1)
                Y(IYBEG) = Y(IYBEG) - A1*A1
  100       CONTINUE
            IYBEG = IYBEG + 2*LENY + 1
            LENY = LENY - 2
            MM = MM - 2
  200   CONTINUE
C       
C       -------------------------------------------------------
C       UPDATE TWO COLUMNS OF Y AT A TIME,  EXCEPT THE DIAGONAL 
C       ELEMENT.
C       NOTE: THE DIAGONAL ELEMENT OF THE ODD COLUMN HAS
C             BEEN COMPUTED, SO WE COMPUTE THE SAME NUMBER OF
C             ELEMENTS FOR THE TWO COLUMNS.
C       -------------------------------------------------------
C
        MM = M
        IYBEG = 1
        LENY = LDY - 1 
C
        DO  600  J = 1, QQ-1, 2
C
            IYBEG1 = IYBEG 
            IYBEG2 = IYBEG + LENY
C
            DO  400  K = 1, N-1, 2
C
C               ---------------------------------
C               TWO COLUMNS UPDATING TWO COLUMNS.
C               ---------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                A1 = X(I1)
                A2 = X(I2)
                A9  = X(I1+1)
                A10 = X(I2+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10
C
                DO  300  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B2 * A10
                    Y(IYBEG2+I) = Y2
  300           CONTINUE
C
  400       CONTINUE
C
C           -----------------------------
C           BOUNDARY CODE FOR THE K LOOP.
C           -----------------------------
C
            IF  ( K .EQ. N )  THEN
C
C               --------------------------------
C               ONE COLUMN UPDATING TWO COLUMNS.
C               --------------------------------
C
                I1 = XPNT(K+1) - MM
                A1 = X(I1)
                A9  = X(I1+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9
C
                DO  500  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B1 * A9
                    Y(IYBEG2+I) = Y2
  500           CONTINUE
C
            ENDIF
C
C           -----------------------------------------------
C           PREPARE FOR NEXT PAIR OF COLUMNS TO BE UPDATED.
C           -----------------------------------------------
C
            MM = MM - 2
            IYBEG = IYBEG2 + LENY + 1
            LENY = LENY - 2
C
  600   CONTINUE
C
C       ------------------------------------------------------
C       BOUNDARY CODE FOR J LOOP:  EXECUTED WHENEVER Q IS ODD.
C       ------------------------------------------------------
C
        IF  ( J .EQ. QQ )  THEN
            CALL  SMXPY2  ( MM, N, Y(IYBEG), XPNT, X )
        ENDIF
C
        RETURN
        END
