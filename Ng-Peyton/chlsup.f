C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  December 27, 1994
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratoy
C
C***********************************************************************
C***********************************************************************
C******     CHLSUP .... DENSE CHOLESKY WITHIN SUPERNODE   **************
C***********************************************************************
C***********************************************************************
C
C     PURPOSE - THIS ROUTINE PERFORMS CHOLESKY
C               FACTORIZATION ON THE COLUMNS OF A SUPERNODE
C               THAT HAVE RECEIVED ALL UPDATES FROM COLUMNS
C               EXTERNAL TO THE SUPERNODE.
C
C     INPUT PARAMETERS -
C        M      - NUMBER OF ROWS (LENGTH OF THE FIRST COLUMN).
C        N      - NUMBER OF COLUMNS IN THE SUPERNODE.
C        XPNT   - XPNT(J+1) POINTS ONE LOCATION BEYOND THE END
C                 OF THE J-TH COLUMN OF THE SUPERNODE.
C        X(*)   - CONTAINS THE COLUMNS OF OF THE SUPERNODE TO
C                 BE FACTORED.
C        SMXPY  - EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY.
C
C     OUTPUT PARAMETERS -
C        X(*)   - ON OUTPUT, CONTAINS THE FACTORED COLUMNS OF
C                 THE SUPERNODE.
C        IFLAG  - UNCHANGED IF THERE IS NO ERROR.
C                 =1 IF NONPOSITIVE DIAGONAL ENTRY IS ENCOUNTERED.
C
C***********************************************************************
C
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
      SUBROUTINE  CHLSUP  ( M, N, SPLIT, XPNT, X, MXDIAG, NTINY, 
     &                      IFLAG, MMPYN, SMXPY )
C
C***********************************************************************
C
C     -----------
C     PARAMETERS.
C     -----------
C
      EXTERNAL            MMPYN, SMXPY
C
      INTEGER             M, N, IFLAG
C
      INTEGER             XPNT(*), SPLIT(*)
C
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
      DOUBLE PRECISION    X(*), MXDIAG
      INTEGER             NTINY
C
C     ----------------
C     LOCAL VARIABLES.
C     ----------------
C
      INTEGER             FSTCOL, JBLK  , JPNT  , MM    , NN    ,
     &                    NXTCOL, Q
C
C***********************************************************************
C
        JBLK = 0
        FSTCOL = 1
        MM = M
        JPNT = XPNT(FSTCOL)
C
C       ----------------------------------------
C       FOR EACH BLOCK JBLK IN THE SUPERNODE ...
C       ----------------------------------------
  100   CONTINUE
        IF  ( FSTCOL .LE. N )  THEN
            JBLK = JBLK + 1
            NN = SPLIT(JBLK)
C           ------------------------------------------
C           ... PERFORM PARTIAL CHOLESKY FACTORIZATION
C               ON THE BLOCK.
C           ------------------------------------------
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
            CALL PCHOL ( MM, NN, XPNT(FSTCOL), X, MXDIAG, NTINY,
     &                                         IFLAG, SMXPY )
            IF  ( IFLAG .EQ. 1 )  RETURN
C           ----------------------------------------------
C           ... APPLY THE COLUMNS IN JBLK TO ANY COLUMNS
C               OF THE SUPERNODE REMAINING TO BE COMPUTED.
C           ----------------------------------------------
            NXTCOL = FSTCOL + NN
            Q = N - NXTCOL + 1
            MM = MM - NN
            JPNT = XPNT(NXTCOL)
            IF  ( Q .GT. 0 )  THEN
                CALL  MMPYN ( MM, NN, Q, XPNT(FSTCOL), X, X(JPNT), MM )
            ENDIF
            FSTCOL = NXTCOL
            GO TO 100
        ENDIF
C
        RETURN
        END
