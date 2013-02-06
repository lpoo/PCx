C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  May 26, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C****     FNSPLT ..... COMPUTE FINE PARTITIONING OF SUPERNODES     *****
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE DETERMINES A FINE PARTITIONING OF SUPERNODES
C       WHEN THERE IS A CACHE AVAILABLE ON THE MACHINE.  THE FINE
C       PARTITIONING IS CHOSEN SO THAT DATA RE-USE IS MAXIMIZED.
C       
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
C                           THE SUPERNODE PARTITIONING.
C       XLINDX          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
C                           POINTERS IN THE SUPERNODE INDICES.
C       CACHSZ          -   CACHE SIZE IN KILO BYTES.
C                           IF THERE IS NO CACHE, SET CACHSZ = 0.
C
C   OUTPUT PARAMETERS:
C       SPLIT           -   INTEGER ARRAY OF SIZE NEQNS CONTAINING THE
C                           FINE PARTITIONING.
C
C***********************************************************************
C
        SUBROUTINE  FNSPLT ( NEQNS , NSUPER, XSUPER, XLINDX,
     &                       CACHSZ, SPLIT )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER         CACHSZ, NEQNS , NSUPER
        INTEGER         XSUPER(*), SPLIT(*)
        INTEGER         XLINDX(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER         CACHE , CURCOL, FSTCOL, HEIGHT, KCOL  , 
     1                  KSUP  , LSTCOL, NCOLS , NXTBLK, USED  , 
     1                  WIDTH
C
C *******************************************************************
C
C       --------------------------------------------
C       COMPUTE THE NUMBER OF 8-BYTE WORDS IN CACHE.
C       --------------------------------------------
        IF  ( CACHSZ .LE. 0 )  THEN
            CACHE = 2 000 000 000
        ELSE
            CACHE = ( FLOAT(CACHSZ) * 1024. / 8. ) * 0.9
        ENDIF
C
C       ---------------
C       INITIALIZATION.
C       ---------------
        DO  100  KCOL = 1, NEQNS
            SPLIT(KCOL) = 0
  100   CONTINUE
C
C       ---------------------------
C       FOR EACH SUPERNODE KSUP ...
C       ---------------------------
        DO  1000  KSUP = 1, NSUPER
C           -----------------------
C           ... GET SUPERNODE INFO.
C           -----------------------
            HEIGHT = XLINDX(KSUP+1) - XLINDX(KSUP)
            FSTCOL = XSUPER(KSUP)
            LSTCOL = XSUPER(KSUP+1) - 1
            WIDTH = LSTCOL - FSTCOL + 1
            NXTBLK = FSTCOL
C           --------------------------------------
C           ... UNTIL ALL COLUMNS OF THE SUPERNODE 
C               HAVE BEEN PROCESSED ...
C           --------------------------------------
            CURCOL = FSTCOL - 1
  200       CONTINUE
C               -------------------------------------------
C               ... PLACE THE FIRST COLUMN(S) IN THE CACHE.
C               -------------------------------------------
                CURCOL = CURCOL + 1
                IF  ( CURCOL .LT. LSTCOL )  THEN
                    CURCOL = CURCOL + 1
                    NCOLS = 2
                    USED = 4 * HEIGHT - 1
                    HEIGHT = HEIGHT - 2
                ELSE
                    NCOLS = 1
                    USED = 3 * HEIGHT
                    HEIGHT = HEIGHT - 1
                ENDIF
C
C               --------------------------------------
C               ... WHILE THE CACHE IS NOT FILLED AND
C                   THERE ARE COLUMNS OF THE SUPERNODE 
C                   REMAINING TO BE PROCESSED ...
C               --------------------------------------
  300           CONTINUE
                IF  ( USED+HEIGHT .LT. CACHE  .AND.
     &                CURCOL      .LT. LSTCOL       )  THEN
C                   --------------------------------
C                   ... ADD ANOTHER COLUMN TO CACHE.
C                   --------------------------------
                    CURCOL = CURCOL + 1
                    NCOLS = NCOLS + 1
                    USED = USED + HEIGHT
                    HEIGHT = HEIGHT - 1
                    GO TO 300
                ENDIF
C               -------------------------------------
C               ... RECORD THE NUMBER OF COLUMNS THAT 
C                   FILLED THE CACHE.
C               -------------------------------------
                SPLIT(NXTBLK) = NCOLS
                NXTBLK = NXTBLK + 1
C               --------------------------
C               ... GO PROCESS NEXT BLOCK.
C               --------------------------
                IF  ( CURCOL .LT. LSTCOL )  GO TO 200
 1000   CONTINUE
C
        RETURN
        END
