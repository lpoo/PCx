C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  February 13, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C*************     SYMFC2 ..... SYMBOLIC FACTORIZATION    **************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE: 
C       THIS ROUTINE PERFORMS SUPERNODAL SYMBOLIC FACTORIZATION ON A 
C       REORDERED LINEAR SYSTEM.  IT ASSUMES ACCESS TO THE COLUMNS 
C       COUNTS, SUPERNODE PARTITION, AND SUPERNODAL ELIMINATION TREE
C       ASSOCIATED WITH THE FACTOR MATRIX L.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS
C       (I) ADJLEN      -   LENGTH OF THE ADJACENCY LIST.
C       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1 CONTAINING POINTERS
C                           TO THE ADJACENCY STRUCTURE.
C       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1 CONTAINING
C                           THE ADJACENCY STRUCTURE.
C       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           POSTORDERING.
C       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS CONTAINING THE
C                           INVERSE OF THE POSTORDERING.
C       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) NSUPER      -   NUMBER OF SUPERNODES.
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           FIRST COLUMN OF EACH SUPERNODE.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS TO BE STORED IN
C                           LINDX(*).
C
C   OUTPUT PARAMETERS:
C       (I) XLINDX      -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS 
C                           INTO THE SUBSCRIPT VECTOR.
C       (I) LINDX       -   ARRAY OF LENGTH MAXSUB, CONTAINING THE
C                           COMPRESSED SUBSCRIPTS.
C       (I) XLNZ        -   COLUMN POINTERS FOR L.
C       (I) FLAG        -   ERROR FLAG:
C                               0 - NO ERROR.
C                               1 - INCONSISTANCY IN THE INPUT.
C       
C   WORKING PARAMETERS:
C       (I) MRGLNK      -   ARRAY OF LENGTH NSUPER, CONTAINING THE 
C                           CHILDREN OF EACH SUPERNODE AS A LINKED LIST.
C       (I) RCHLNK      -   ARRAY OF LENGTH NEQNS+1, CONTAINING THE 
C                           CURRENT LINKED LIST OF MERGED INDICES (THE 
C                           "REACH" SET).
C       (I) MARKER      -   ARRAY OF LENGTH NEQNS USED TO MARK INDICES
C                           AS THEY ARE INTRODUCED INTO EACH SUPERNODE'S
C                           INDEX SET.
C
C***********************************************************************
C
      SUBROUTINE  SYMFC2 (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  , 
     &                      INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                      NOFSUB, XLINDX, LINDX , XLNZ  , MRGLNK,
     &                      RCHLNK, MARKER, FLAG    )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             ADJLEN, FLAG  , NEQNS , NOFSUB, NSUPER
        INTEGER             ADJNCY(ADJLEN), COLCNT(NEQNS) ,
     &                      INVP(NEQNS)   , MARKER(NEQNS) ,
     &                      MRGLNK(NSUPER), LINDX(NOFSUB) , 
     &                      PERM(NEQNS)   , RCHLNK(0:NEQNS), 
     &                      SNODE(NEQNS)  , XSUPER(NSUPER+1)
        INTEGER             XADJ(NEQNS+1) , XLINDX(NSUPER+1),
     &                      XLNZ(NEQNS+1)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             FSTCOL, HEAD  , I     , JNZBEG, JNZEND,
     &                      JPTR  , JSUP  , JWIDTH, KNZ   , KNZBEG,
     &                      KNZEND, KPTR  , KSUP  , LENGTH, LSTCOL,
     &                      NEWI  , NEXTI , NODE  , NZBEG , NZEND ,
     &                      PCOL  , PSUP  , POINT , TAIL  , WIDTH
C
C***********************************************************************
C
        FLAG = 0
        IF  ( NEQNS .LE. 0 )  RETURN
C
C       ---------------------------------------------------
C       INITIALIZATIONS ...
C           NZEND  : POINTS TO THE LAST USED SLOT IN LINDX.
C           TAIL   : END OF LIST INDICATOR 
C                    (IN RCHLNK(*), NOT MRGLNK(*)).
C           MRGLNK : CREATE EMPTY LISTS.
C           MARKER : "UNMARK" THE INDICES.
C       ---------------------------------------------------
        NZEND = 0
        HEAD = 0
        TAIL = NEQNS + 1
        POINT = 1
        DO  50  I = 1, NEQNS
            MARKER(I) = 0
            XLNZ(I) = POINT
            POINT = POINT + COLCNT(I)
   50   CONTINUE
        XLNZ(NEQNS+1) = POINT
        POINT = 1
        DO  100  KSUP = 1, NSUPER
            MRGLNK(KSUP) = 0
            FSTCOL = XSUPER(KSUP)
            XLINDX(KSUP) = POINT
            POINT = POINT + COLCNT(FSTCOL)
  100   CONTINUE
        XLINDX(NSUPER+1) = POINT
C
C       ---------------------------
C       FOR EACH SUPERNODE KSUP ... 
C       ---------------------------
        DO  1000  KSUP = 1, NSUPER
C
C           ---------------------------------------------------------
C           INITIALIZATIONS ...
C               FSTCOL : FIRST COLUMN OF SUPERNODE KSUP.
C               LSTCOL : LAST COLUMN OF SUPERNODE KSUP.
C               KNZ    : WILL COUNT THE NONZEROS OF L IN COLUMN KCOL.
C               RCHLNK : INITIALIZE EMPTY INDEX LIST FOR KCOL.
C           ---------------------------------------------------------
            FSTCOL = XSUPER(KSUP)
            LSTCOL = XSUPER(KSUP+1) - 1
            WIDTH  = LSTCOL - FSTCOL + 1
            LENGTH = COLCNT(FSTCOL)
            KNZ = 0
            RCHLNK(HEAD) = TAIL
            JSUP = MRGLNK(KSUP)
C
C           -------------------------------------------------
C           IF KSUP HAS CHILDREN IN THE SUPERNODAL E-TREE ...
C           -------------------------------------------------
            IF  ( JSUP .GT. 0 )  THEN
C               ---------------------------------------------
C               COPY THE INDICES OF THE FIRST CHILD JSUP INTO 
C               THE LINKED LIST, AND MARK EACH WITH THE VALUE 
C               KSUP.
C               ---------------------------------------------
                JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
                JNZBEG = XLINDX(JSUP) + JWIDTH
                JNZEND = XLINDX(JSUP+1) - 1
                DO  200  JPTR = JNZEND, JNZBEG, -1
                    NEWI = LINDX(JPTR)
                    KNZ = KNZ+1
                    MARKER(NEWI) = KSUP
                    RCHLNK(NEWI) = RCHLNK(HEAD)
                    RCHLNK(HEAD) = NEWI
  200           CONTINUE
C               ------------------------------------------
C               FOR EACH SUBSEQUENT CHILD JSUP OF KSUP ...
C               ------------------------------------------
                JSUP = MRGLNK(JSUP)
  300           CONTINUE
                IF  ( JSUP .NE. 0  .AND.  KNZ .LT. LENGTH )  THEN
C                   ----------------------------------------
C                   MERGE THE INDICES OF JSUP INTO THE LIST,
C                   AND MARK NEW INDICES WITH VALUE KSUP.
C                   ----------------------------------------
                    JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
                    JNZBEG = XLINDX(JSUP) + JWIDTH
                    JNZEND = XLINDX(JSUP+1) - 1
                    NEXTI = HEAD
                    DO  500  JPTR = JNZBEG, JNZEND
                        NEWI = LINDX(JPTR)
  400                   CONTINUE
                            I = NEXTI
                            NEXTI = RCHLNK(I)
                            IF  ( NEWI .GT. NEXTI )  GO TO 400
                        IF  ( NEWI .LT. NEXTI )  THEN
                            KNZ = KNZ+1
                            RCHLNK(I) = NEWI
                            RCHLNK(NEWI) = NEXTI
                            MARKER(NEWI) = KSUP
                            NEXTI = NEWI
                        ENDIF
  500               CONTINUE
                    JSUP = MRGLNK(JSUP)
                    GO TO 300
                ENDIF
            ENDIF
C           ---------------------------------------------------
C           STRUCTURE OF A(*,FSTCOL) HAS NOT BEEN EXAMINED YET.  
C           "SORT" ITS STRUCTURE INTO THE LINKED LIST,
C           INSERTING ONLY THOSE INDICES NOT ALREADY IN THE
C           LIST.
C           ---------------------------------------------------
            IF  ( KNZ .LT. LENGTH )  THEN
                NODE = PERM(FSTCOL)
                KNZBEG = XADJ(NODE)
                KNZEND = XADJ(NODE+1) - 1
                DO  700  KPTR = KNZBEG, KNZEND
                    NEWI = ADJNCY(KPTR)
                    NEWI = INVP(NEWI)
                    IF  ( NEWI .GT. FSTCOL  .AND.
     &                    MARKER(NEWI) .NE. KSUP )  THEN
C                       --------------------------------
C                       POSITION AND INSERT NEWI IN LIST
C                       AND MARK IT WITH KCOL.
C                       --------------------------------
                        NEXTI = HEAD
  600                   CONTINUE
                            I = NEXTI
                            NEXTI = RCHLNK(I)
                            IF  ( NEWI .GT. NEXTI )  GO TO 600
                        KNZ = KNZ + 1
                        RCHLNK(I) = NEWI
                        RCHLNK(NEWI) = NEXTI
                        MARKER(NEWI) = KSUP
                    ENDIF
  700           CONTINUE
            ENDIF
C           ------------------------------------------------------------
C           IF KSUP HAS NO CHILDREN, INSERT FSTCOL INTO THE LINKED LIST.
C           ------------------------------------------------------------
            IF  ( RCHLNK(HEAD) .NE. FSTCOL )  THEN
                RCHLNK(FSTCOL) = RCHLNK(HEAD)
                RCHLNK(HEAD) = FSTCOL
                KNZ = KNZ + 1
            ENDIF
C
C           --------------------------------------------
C           COPY INDICES FROM LINKED LIST INTO LINDX(*).
C           --------------------------------------------
            NZBEG = NZEND + 1
            NZEND = NZEND + KNZ
            IF  ( NZEND+1 .NE. XLINDX(KSUP+1) )  GO TO 8000
            I = HEAD
            DO  800  KPTR = NZBEG, NZEND
                I = RCHLNK(I)
                LINDX(KPTR) = I
  800       CONTINUE
C
C           ---------------------------------------------------
C           IF KSUP HAS A PARENT, INSERT KSUP INTO ITS PARENT'S 
C           "MERGE" LIST.
C           ---------------------------------------------------
            IF  ( LENGTH .GT. WIDTH )  THEN
                PCOL = LINDX ( XLINDX(KSUP) + WIDTH )
                PSUP = SNODE(PCOL)
                MRGLNK(KSUP) = MRGLNK(PSUP)
                MRGLNK(PSUP) = KSUP
            ENDIF
C
 1000   CONTINUE
C
        RETURN
C
C       -----------------------------------------------
C       INCONSISTENCY IN DATA STRUCTURE WAS DISCOVERED.
C       -----------------------------------------------
 8000   CONTINUE
        FLAG = -2
        RETURN
C
      END
