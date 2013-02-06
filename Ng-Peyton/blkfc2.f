C***********************************************************************
C***********************************************************************
C
C   Version:        0.3
C   Last modified:  March 6, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratoy
C
C***********************************************************************
C***********************************************************************
C*********     BLKFC2 .....  BLOCK GENERAL SPARSE CHOLESKY     *********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE FACTORS A SPARSE POSITIVE DEFINITE MATRIX.
C       THE COMPUTATION IS ORGANIZED AROUND KERNELS THAT PERFORM
C       SUPERNODE-TO-SUPERNODE UPDATES, I.E., BLOCK-TO-BLOCK UPDATES.
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
C       TMPSIZ          -   SIZE OF TEMPORARY WORKING STORAGE.
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
C
C   WORKING PARAMETERS:
C       LINK            -   LINKS TOGETHER THE SUPERNODES IN A SUPERNODE
C                           ROW.
C       LENGTH          -   LENGTH OF THE ACTIVE PORTION OF EACH 
C                           SUPERNODE.
C       INDMAP          -   VECTOR OF SIZE NEQNS INTO WHICH THE GLOBAL
C                           INDICES ARE SCATTERED.
C       RELIND          -   MAPS LOCATIONS IN THE UPDATING COLUMNS TO 
C                           THE CORRESPONDING LOCATIONS IN THE UPDATED 
C                           COLUMNS.  (RELIND IS GATHERED FROM INDMAP).
C       TEMP            -   REAL VECTOR FOR ACCUMULATING UPDATES.  MUST
C                           ACCOMODATE ALL COLUMNS OF A SUPERNODE. 
C       
C***********************************************************************
C
      SUBROUTINE  BLKFC2 (  NSUPER, XSUPER, SNODE , SPLIT , XLINDX,
     &                      LINDX , XLNZ  , LNZ   , LINK  , LENGTH,
     &                      INDMAP, RELIND, TMPSIZ, TEMP  , IFLAG ,
     &                      MMPYN , SMXPY                           )
C
C*********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        EXTERNAL            MMPYN , SMXPY
        INTEGER             XLINDX(*)     , XLNZ(*)
        INTEGER             INDMAP(*)     , LENGTH(*)     ,
     &                      LINDX(*)      , LINK(*)       ,
     &                      RELIND(*)     , SNODE(*)      ,
     &                      SPLIT(*)      , XSUPER(*)
        INTEGER             IFLAG , NSUPER, TMPSIZ
        DOUBLE PRECISION    LNZ(*)        , TEMP(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
C
        INTEGER             FJCOL , FKCOL , I     , ILEN  , ILPNT ,
     &                      INDDIF, JLEN  , JLPNT , JSUP  , JXPNT ,
     &                      KFIRST, KLAST , KLEN  , KLPNT , KSUP  ,
     &                      KXPNT , LJCOL , NCOLUP, NJCOLS, NKCOLS,
     &                      NXKSUP, NXTCOL, NXTSUP, STORE

CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
        DOUBLE PRECISION MXDIAG
        INTEGER NTINY
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
C
C*********************************************************************
C
        IFLAG = 0
        NTINY = 0
C
C       -----------------------------------------------------------
C       INITIALIZE EMPTY ROW LISTS IN LINK(*) AND ZERO OUT TEMP(*).
C       -----------------------------------------------------------
        DO  100  JSUP = 1, NSUPER
            LINK(JSUP) = 0
  100   CONTINUE
        DO  200  I = 1, TMPSIZ
            TEMP(I) = 0.0D+00
  200   CONTINUE
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
C       COMPUTE MAXIMUM DIAGONAL ELEMENT IN INPUT MATRIX
        MXDIAG = 0.D0
        DO 201 I = 1, XSUPER(NSUPER+1)-1
          FJCOL = XLNZ(I)
          MXDIAG = MAX(MXDIAG, LNZ(FJCOL))
 201    CONTINUE
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
C
C       ---------------------------
C       FOR EACH SUPERNODE JSUP ...
C       ---------------------------
        DO  600  JSUP = 1, NSUPER
C
C           ------------------------------------------------
C           FJCOL  ...  FIRST COLUMN OF SUPERNODE JSUP.
C           LJCOL  ...  LAST COLUMN OF SUPERNODE JSUP.
C           NJCOLS ...  NUMBER OF COLUMNS IN SUPERNODE JSUP.
C           JLEN   ...  LENGTH OF COLUMN FJCOL.
C           JXPNT  ...  POINTER TO INDEX OF FIRST
C                       NONZERO IN COLUMN FJCOL.
C           ------------------------------------------------
            FJCOL  = XSUPER(JSUP)
            NJCOLS = XSUPER(JSUP+1) - FJCOL
            LJCOL  = FJCOL + NJCOLS - 1
            JLEN   = XLNZ(FJCOL+1) - XLNZ(FJCOL)
            JXPNT  = XLINDX(JSUP)

C            print *, 'Super Node: ', JSUP, ' first: ', FJCOL, 
C     .           ' last: ', LJCOL

            
C
C
C           -----------------------------------------------------
C           SET UP INDMAP(*) TO MAP THE ENTRIES IN UPDATE COLUMNS
C           TO THEIR CORRESPONDING POSITIONS IN UPDATED COLUMNS, 
C           RELATIVE THE THE BOTTOM OF EACH UPDATED COLUMN.
C           -----------------------------------------------------
            CALL  LDINDX ( JLEN, LINDX(JXPNT), INDMAP )
C
C           -----------------------------------------
C           FOR EVERY SUPERNODE KSUP IN ROW(JSUP) ...
C           -----------------------------------------
            KSUP = LINK(JSUP)
  300       IF  ( KSUP .GT. 0 )  THEN
                NXKSUP = LINK(KSUP)
C
C               -------------------------------------------------------
C               GET INFO ABOUT THE CMOD(JSUP,KSUP) UPDATE.
C
C               FKCOL  ...  FIRST COLUMN OF SUPERNODE KSUP.
C               NKCOLS ...  NUMBER OF COLUMNS IN SUPERNODE KSUP.
C               KLEN   ...  LENGTH OF ACTIVE PORTION OF COLUMN FKCOL.
C               KXPNT  ...  POINTER TO INDEX OF FIRST NONZERO IN ACTIVE
C                           PORTION OF COLUMN FJCOL.
C               -------------------------------------------------------
                FKCOL = XSUPER(KSUP)
                NKCOLS = XSUPER(KSUP+1) - FKCOL
                KLEN = LENGTH(KSUP)
                KXPNT = XLINDX(KSUP+1) - KLEN
C
C               -------------------------------------------
C               PERFORM CMOD(JSUP,KSUP), WITH SPECIAL CASES
C               HANDLED DIFFERENTLY.
C               -------------------------------------------
C
                IF  ( KLEN .NE. JLEN )  THEN
C
C                   -------------------------------------------
C                   SPARSE CMOD(JSUP,KSUP).
C
C                   NCOLUP ... NUMBER OF COLUMNS TO BE UPDATED.
C                   -------------------------------------------
C
                    DO  400  I = 0, KLEN-1
                        NXTCOL = LINDX(KXPNT+I)
                        IF  ( NXTCOL .GT. LJCOL )  GO TO 500
  400               CONTINUE
                    I = KLEN
  500               CONTINUE
                    NCOLUP = I
C
                    IF  ( NKCOLS .EQ. 1 )  THEN
C
C                       ----------------------------------------------
C                       UPDATING TARGET SUPERNODE BY TRIVIAL
C                       SUPERNODE (WITH ONE COLUMN).
C
C                       KLPNT  ...  POINTER TO FIRST NONZERO IN ACTIVE
C                                   PORTION OF COLUMN FKCOL.
C                       ----------------------------------------------
                        KLPNT = XLNZ(FKCOL+1) - KLEN
                        CALL  MMPYI ( KLEN, NCOLUP, LINDX(KXPNT),
     &                                LNZ(KLPNT), XLNZ, LNZ, INDMAP )
C
                    ELSE
C
C                       --------------------------------------------
C                       KFIRST ...  FIRST INDEX OF ACTIVE PORTION OF
C                                   SUPERNODE KSUP (FIRST COLUMN TO
C                                   BE UPDATED).
C                       KLAST  ...  LAST INDEX OF ACTIVE PORTION OF
C                                   SUPERNODE KSUP.
C                       --------------------------------------------
C
                        KFIRST = LINDX(KXPNT)
                        KLAST  = LINDX(KXPNT+KLEN-1)
                        INDDIF = INDMAP(KFIRST) - INDMAP(KLAST)
C
                        IF  ( INDDIF .LT. KLEN )  THEN
C
C                           ---------------------------------------
C                           DENSE CMOD(JSUP,KSUP).
C
C                           ILPNT  ...  POINTER TO FIRST NONZERO IN
C                                       COLUMN KFIRST.
C                           ILEN   ...  LENGTH OF COLUMN KFIRST.
C                           ---------------------------------------
                            ILPNT = XLNZ(KFIRST)
                            ILEN = XLNZ(KFIRST+1) - ILPNT
                            CALL  MMPY ( KLEN, NKCOLS, NCOLUP,
     &                                   SPLIT(FKCOL), XLNZ(FKCOL),
     &                                   LNZ, LNZ(ILPNT), ILEN, MMPYN  )
C
                        ELSE
C
C                           -------------------------------
C                           GENERAL SPARSE CMOD(JSUP,KSUP).
C                           COMPUTE CMOD(JSUP,KSUP) UPDATE
C                           IN WORK STORAGE.
C                           -------------------------------
                            STORE = KLEN * NCOLUP - NCOLUP * 
     &                              (NCOLUP-1) / 2
                            IF  ( STORE .GT. TMPSIZ )  THEN
                                IFLAG = -2
                                RETURN
                            ENDIF
                            CALL  MMPY ( KLEN, NKCOLS, NCOLUP,
     &                                   SPLIT(FKCOL), XLNZ(FKCOL),
     &                                   LNZ, TEMP, KLEN, MMPYN  )
C                           ----------------------------------------
C                           GATHER INDICES OF KSUP RELATIVE TO JSUP.
C                           ----------------------------------------
                            CALL  IGATHR ( KLEN, LINDX(KXPNT),
     &                                     INDMAP, RELIND )
C                           --------------------------------------
C                           INCORPORATE THE CMOD(JSUP,KSUP) BLOCK
C                           UPDATE INTO THE TO APPROPRIATE COLUMNS
C                           OF L.
C                           --------------------------------------
                            CALL  ASSMB ( KLEN, NCOLUP, TEMP, RELIND,
     &                                    XLNZ(FJCOL), LNZ, JLEN )
C
                        ENDIF
C
                    ENDIF
C
                ELSE
C
C                   ----------------------------------------------
C                   DENSE CMOD(JSUP,KSUP).
C                   JSUP AND KSUP HAVE IDENTICAL STRUCTURE.
C
C                   JLPNT  ...  POINTER TO FIRST NONZERO IN COLUMN
C                               FJCOL.
C                   ----------------------------------------------
                    JLPNT = XLNZ(FJCOL)
                    CALL  MMPY ( KLEN, NKCOLS, NJCOLS, SPLIT(FKCOL),
     &                           XLNZ(FKCOL), LNZ, LNZ(JLPNT), JLEN,
     &                           MMPYN )
                    NCOLUP = NJCOLS
                    IF  ( KLEN .GT. NJCOLS )  THEN
                        NXTCOL = LINDX(JXPNT+NJCOLS)
                    ENDIF
C
                ENDIF
C
C               ------------------------------------------------
C               LINK KSUP INTO LINKED LIST OF THE NEXT SUPERNODE
C               IT WILL UPDATE AND DECREMENT KSUP'S ACTIVE
C               LENGTH.
C               ------------------------------------------------
                IF  ( KLEN .GT. NCOLUP )  THEN
                    NXTSUP = SNODE(NXTCOL)
                    LINK(KSUP) = LINK(NXTSUP)
                    LINK(NXTSUP) = KSUP
                    LENGTH(KSUP) = KLEN - NCOLUP
                ELSE
                    LENGTH(KSUP) = 0
                ENDIF
C
C               -------------------------------
C               NEXT UPDATING SUPERNODE (KSUP).
C               -------------------------------
                KSUP = NXKSUP
                GO TO 300
C
            ENDIF
C
C           ----------------------------------------------
C           APPLY PARTIAL CHOLESKY TO THE COLUMNS OF JSUP.
C           ----------------------------------------------
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
            CALL CHLSUP ( JLEN, NJCOLS, SPLIT(FJCOL), XLNZ(FJCOL), LNZ,
     &                    MXDIAG, NTINY, IFLAG, MMPYN, SMXPY )
            IF  ( IFLAG .NE. 0 )  THEN
                IFLAG = -1
                RETURN
            ENDIF
C
C           -----------------------------------------------
C           INSERT JSUP INTO LINKED LIST OF FIRST SUPERNODE
C           IT WILL UPDATE.
C           -----------------------------------------------
            IF  ( JLEN .GT. NJCOLS )  THEN
                NXTCOL = LINDX(JXPNT+NJCOLS)
                NXTSUP = SNODE(NXTCOL)
                LINK(JSUP) = LINK(NXTSUP)
                LINK(NXTSUP) = JSUP
                LENGTH(JSUP) = JLEN - NJCOLS
            ELSE
                LENGTH(JSUP) = 0
            ENDIF
C
  600   CONTINUE
C
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
        IF(NTINY .NE. 0) WRITE(6,699) NTINY
 699    FORMAT(1X,' FOUND ',I6,' TINY DIAGONALS; REPLACED WITH INF')
C
C SET IFLAG TO -1 TO INDICATE PRESENCE OF TINY DIAGONALS
C
	IF(NTINY .NE. 0) IFLAG = -1
CxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPCxPC
        RETURN
      END
