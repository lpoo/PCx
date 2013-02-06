C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  March 6, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C*********     BLKFCT .....  BLOCK GENERAL SPARSE CHOLESKY     *********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE CALLS THE BLOCK GENERAL SPARSE CHOLESKY ROUTINE,
C       BLKFC2.
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       SNODE           -   MAPS EACH COLUMN TO THE SUPERNODE CONTAINING
C                           IT.
C       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT
C                           INTO CACHE.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE (INCLUDING
C                           THE DIAGONAL ELEMENTS).
C       (XLNZ,LNZ)      -   ON INPUT, CONTAINS MATRIX TO BE FACTORED.
C       IWSIZ           -   SIZE OF INTEGER WORKING STORAGE
C       TMPSIZ          -   SIZE OF FLOATING POINT WORKING STORAGE.
C       MMPYN           -   EXTERNAL ROUTINE: MATRIX-MATRIX MULTIPLY.
C       SMXPY           -   EXTERNAL ROUTINE: MATRIX-VECTOR MULTIPLY.
C
C   OUTPUT PARAMETERS:
C       LNZ             -   ON OUTPUT, CONTAINS CHOLESKY FACTOR.
C       IFLAG           -   ERROR FLAG.
C                               0: SUCCESSFUL FACTORIZATION.
C                              -1: NONPOSITIVE DIAGONAL ENCOUNTERED,
C                                  MATRIX IS NOT POSITIVE DEFINITE.
C                              -2: INSUFFICIENT WORKING STORAGE 
C                                  [TEMP(*)].
C                              -3: INSUFFICIENT WORKING STORAGE 
C                                  [IWORK(*)].
C
C   WORKING PARAMETERS:
C       IWORK           -   INTEGER WORKING STORAGE OF LENGTH 
C                           2*NEQNS + 2*NSUPER.
C       TMPVEC          -   DOUBLE PRECISION WORKING STORAGE OF LENGTH
C                           NEQNS.
C       
C***********************************************************************
C
      SUBROUTINE  BLKFCT (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPYN , 
     &                      SMXPY                                   )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        EXTERNAL            MMPYN , SMXPY
        INTEGER             XLINDX(*)     , XLNZ(*)
        INTEGER             IWORK(*)      , LINDX(*)      , 
     &                      SNODE(*)      , SPLIT(*)      , 
     &                      XSUPER(*)
        INTEGER             IFLAG , IWSIZ , NEQNS , NSUPER, TMPSIZ
        DOUBLE PRECISION    LNZ(*)        , TMPVEC(*)
C
C*********************************************************************
C
        IFLAG = 0
        IF  ( IWSIZ .LT. 2*NEQNS+2*NSUPER )  THEN
            IFLAG = -3
            RETURN
        ENDIF
        CALL  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , XLINDX,
     &                  LINDX , XLNZ  , LNZ   , 
     &                  IWORK(1)                      ,
     &                  IWORK(NSUPER+1)               ,
     &                  IWORK(2*NSUPER+1)             ,
     &                  IWORK(2*NSUPER+NEQNS+1)       ,
     &                  TMPSIZ, TMPVEC, IFLAG , MMPYN , SMXPY   )
        RETURN
      END
