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
C*************     MMPY8  .... MATRIX-MATRIX MULTIPLY     **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE -
C       THIS ROUTINE PERFORMS A MATRIX-MATRIX MULTIPLY, Y = Y + XA,
C       ASSUMING DATA STRUCTURES USED IN SOME OF OUR SPARSE CHOLESKY
C       CODES.
C
C       LOOP UNROLLING: LEVEL 8 UPDATING TWO COLUMNS AT A TIME
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
        SUBROUTINE  MMPY8  (  M     , N     , Q     , XPNT  , X     ,
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
        INTEGER               I1    , I2    , I3    , I4    , I5    ,
     &                        I6    , I7    , I8
        INTEGER               IYBEG , IYBEG1, IYBEG2, LENY  , MM    
        DOUBLE PRECISION      A1    , A2    , A3    , A4    , A5    ,
     &                        A6    , A7    , A8    , A9    , A10   ,
     &                        A11   , A12   , A13   , A14   , A15   ,
     &                        A16
        DOUBLE PRECISION      B1    , B2    , B3    , B4    , B5    ,
     &                        B6    , B7    , B8    , Y1    , Y2
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
        DO  3000  J = 1, QQ-1, 2
C
            IYBEG1 = IYBEG 
            IYBEG2 = IYBEG + LENY
C
            DO  400  K = 1, N-7, 8
C
C               -----------------------------------
C               EIGHT COLUMNS UPDATING TWO COLUMNS.
C               -----------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                I4 = XPNT(K+4) - MM
                I5 = XPNT(K+5) - MM
                I6 = XPNT(K+6) - MM
                I7 = XPNT(K+7) - MM
                I8 = XPNT(K+8) - MM
                A1 = X(I1)
                A2 = X(I2)
                A3 = X(I3)
                A4 = X(I4)
                A5 = X(I5)
                A6 = X(I6)
                A7 = X(I7)
                A8 = X(I8)
                A9  = X(I1+1)
                A10 = X(I2+1)
                A11 = X(I3+1)
                A12 = X(I4+1)
                A13 = X(I5+1)
                A14 = X(I6+1)
                A15 = X(I7+1)
                A16 = X(I8+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - 
     &              A6*A14 - A7*A15 - A8*A16
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - 
     &              A14*A14 - A15*A15 - A16*A16
C
                DO  300  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    B4 = X(I4+I)
                    Y2 =  Y2 - B3 * A11
                    Y1 =  Y1 - B4 * A4
                    B5 = X(I5+I)
                    Y2 =  Y2 - B4 * A12
                    Y1 =  Y1 - B5 * A5
                    B6 = X(I6+I)
                    Y2 =  Y2 - B5 * A13
                    Y1 =  Y1 - B6 * A6
                    B7 = X(I7+I)
                    Y2 =  Y2 - B6 * A14
                    Y1 = Y1 - B7 * A7
                    B8 = X(I8+I)
                    Y2 = Y2 - B7 * A15
                    Y1 = Y1 - B8 * A8       
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B8 * A16
                    Y(IYBEG2+I) = Y2               
  300           CONTINUE
C
  400       CONTINUE
C
C           -----------------------------
C           BOUNDARY CODE FOR THE K LOOP.
C           -----------------------------
C
            GO TO ( 2000, 1700, 1500, 1300,
     &              1100,  900,  700,  500  ), N-K+2
C
  500       CONTINUE
C
C               -----------------------------------
C               SEVEN COLUMNS UPDATING TWO COLUMNS.
C               -----------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                I4 = XPNT(K+4) - MM
                I5 = XPNT(K+5) - MM
                I6 = XPNT(K+6) - MM
                I7 = XPNT(K+7) - MM
                A1 = X(I1) 
                A2 = X(I2) 
                A3 = X(I3) 
                A4 = X(I4) 
                A5 = X(I5) 
                A6 = X(I6) 
                A7 = X(I7) 
                A9  = X(I1+1) 
                A10 = X(I2+1)
                A11 = X(I3+1)
                A12 = X(I4+1)
                A13 = X(I5+1)
                A14 = X(I6+1)
                A15 = X(I7+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - 
     &              A6*A14 - A7*A15
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - 
     &              A14*A14 - A15*A15
C
                DO  600  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    B4 = X(I4+I)
                    Y2 =  Y2 - B3 * A11
                    Y1 =  Y1 - B4 * A4
                    B5 = X(I5+I)
                    Y2 =  Y2 - B4 * A12
                    Y1 =  Y1 - B5 * A5
                    B6 = X(I6+I)
                    Y2 =  Y2 - B5 * A13
                    Y1 =  Y1 - B6 * A6
                    B7 = X(I7+I)
                    Y2 =  Y2 - B6 * A14
                    Y1 = Y1 - B7 * A7
                    Y(IYBEG1+I) = Y1
                    Y2 = Y2 - B7 * A15
                    Y(IYBEG2+I) = Y2               
  600           CONTINUE
C
                GO TO 2000
C
  700       CONTINUE
C
C               ---------------------------------
C               SIX COLUMNS UPDATING TWO COLUMNS.
C               ---------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                I4 = XPNT(K+4) - MM
                I5 = XPNT(K+5) - MM
                I6 = XPNT(K+6) - MM
                A1 = X(I1)
                A2 = X(I2)
                A3 = X(I3)
                A4 = X(I4)
                A5 = X(I5)
                A6 = X(I6)
                A9  = X(I1+1)
                A10 = X(I2+1)
                A11 = X(I3+1)
                A12 = X(I4+1)
                A13 = X(I5+1)
                A14 = X(I6+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13 - 
     &              A6*A14
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13 - 
     &              A14*A14
C
                DO  800  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    B4 = X(I4+I)
                    Y2 =  Y2 - B3 * A11
                    Y1 =  Y1 - B4 * A4
                    B5 = X(I5+I)
                    Y2 =  Y2 - B4 * A12
                    Y1 =  Y1 - B5 * A5
                    B6 = X(I6+I)
                    Y2 =  Y2 - B5 * A13
                    Y1 =  Y1 - B6 * A6
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B6 * A14
                    Y(IYBEG2+I) = Y2               
  800           CONTINUE
C
                GO TO 2000
C
  900       CONTINUE
C
C               ----------------------------------
C               FIVE COLUMNS UPDATING TWO COLUMNS.
C               ----------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                I4 = XPNT(K+4) - MM
                I5 = XPNT(K+5) - MM
                A1 = X(I1)
                A2 = X(I2)
                A3 = X(I3)
                A4 = X(I4)
                A5 = X(I5)
                A9  = X(I1+1)
                A10 = X(I2+1)
                A11 = X(I3+1)
                A12 = X(I4+1)
                A13 = X(I5+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11 - A4*A12 - A5*A13
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11 - A12*A12 - A13*A13
C
                DO  1000  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    B4 = X(I4+I)
                    Y2 =  Y2 - B3 * A11
                    Y1 =  Y1 - B4 * A4
                    B5 = X(I5+I)
                    Y2 =  Y2 - B4 * A12
                    Y1 =  Y1 - B5 * A5
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B5 * A13
                    Y(IYBEG2+I) = Y2               
 1000           CONTINUE
C
                GO TO 2000
C
 1100       CONTINUE
C
C               ----------------------------------
C               FOUR COLUMNS UPDATING TWO COLUMNS.
C               ----------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                I4 = XPNT(K+4) - MM
                A1 = X(I1)
                A2 = X(I2)
                A3 = X(I3)
                A4 = X(I4)
                A9  = X(I1+1)
                A10 = X(I2+1)
                A11 = X(I3+1)
                A12 = X(I4+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11 - A4*A12
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11 - A12*A12
C
                DO  1200  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    B4 = X(I4+I)
                    Y2 =  Y2 - B3 * A11
                    Y1 =  Y1 - B4 * A4
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B4 * A12
                    Y(IYBEG2+I) = Y2               
 1200           CONTINUE
C
                GO TO 2000
C
 1300       CONTINUE
C
C               -----------------------------------
C               THREE COLUMNS UPDATING TWO COLUMNS.
C               -----------------------------------
C
                I1 = XPNT(K+1) - MM
                I2 = XPNT(K+2) - MM
                I3 = XPNT(K+3) - MM
                A1 = X(I1)
                A2 = X(I2)
                A3 = X(I3)
                A9  = X(I1+1)
                A10 = X(I2+1)
                A11 = X(I3+1)
C
                Y(IYBEG1+1) =  Y(IYBEG1+1) -
     &              A1*A9 - A2*A10 - A3*A11
C
                Y(IYBEG2+1) =  Y(IYBEG2+1) -
     &              A9*A9 - A10*A10 - A11*A11
C
                DO  1400  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    B2 = X(I2+I)
                    Y2 =  Y2 - B1 * A9
                    Y1 =  Y1 - B2 * A2
                    B3 = X(I3+I)
                    Y2 =  Y2 - B2 * A10
                    Y1 =  Y1 - B3 * A3
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B3 * A11
                    Y(IYBEG2+I) = Y2               
 1400           CONTINUE
C
                GO TO 2000
C
 1500       CONTINUE
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
                DO  1600  I = 2, MM-1
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
 1600           CONTINUE
C
                GO TO 2000
C
 1700       CONTINUE
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
                DO  1800  I = 2, MM-1
                    Y1 = Y(IYBEG1+I)
                    B1 = X(I1+I)
                    Y1 =  Y1 - B1 * A1
                    Y2 = Y(IYBEG2+I)
                    Y(IYBEG1+I) = Y1
                    Y2 =  Y2 - B1 * A9
                    Y(IYBEG2+I) = Y2               
 1800           CONTINUE
C
                GO TO 2000
C
C           -----------------------------------------------
C           PREPARE FOR NEXT PAIR OF COLUMNS TO BE UPDATED.
C           -----------------------------------------------
C
 2000       CONTINUE
            MM = MM - 2
            IYBEG = IYBEG2 + LENY + 1
            LENY = LENY - 2
C
 3000   CONTINUE
C
C       -----------------------------------------------------
C       BOUNDARY CODE FOR J LOOP:  EXECUTED WHENVER Q IS ODD.
C       -----------------------------------------------------
C
        IF  ( J .EQ. QQ )  THEN
            CALL  SMXPY8  ( MM, N, Y(IYBEG), XPNT, X )
        ENDIF
C
        RETURN
        END
