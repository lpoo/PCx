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
C****************    FSUP1 ..... FIND SUPERNODES #1    *****************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE IS THE FIRST OF TWO ROUTINES FOR FINDING A
C       MAXIMAL SUPERNODE PARTITION.  IT RETURNS ONLY THE NUMBER OF
C       SUPERNODES NSUPER AND THE SUPERNODE MEMBERSHIP VECTOR SNODE(*), 
C       WHICH IS OF LENGTH NEQNS.  THE VECTORS OF LENGTH NSUPER ARE 
C       COMPUTED SUBSEQUENTLY BY THE COMPANION ROUTINE FSUP2.
C
C   METHOD AND ASSUMPTIONS:
C       THIS ROUTINE USES THE ELIMINATION TREE AND THE FACTOR COLUMN 
C       COUNTS TO COMPUTE THE SUPERNODE PARTITION; IT ALSO ASSUMES A 
C       POSTORDERING OF THE ELIMINATION TREE.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS.
C       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
C       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           FACTOR COLUMN COUNTS: I.E., THE NUMBER OF 
C                           NONZERO ENTRIES IN EACH COLUMN OF L
C                           (INCLUDING THE DIAGONAL ENTRY).
C
C   OUTPUT PARAMETERS:
C       (I) NOFSUB      -   NUMBER OF SUBSCRIPTS.
C       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
C
C   FIRST CREATED ON    JANUARY 18, 1992.
C   LAST UPDATED ON     NOVEMBER 11, 1994.
C
C***********************************************************************
C
      SUBROUTINE  FSUP1  (  NEQNS , ETPAR , COLCNT, NOFSUB, NSUPER,
     &                      SNODE                                   )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             NEQNS , NOFSUB, NSUPER
        INTEGER             COLCNT(*)     , ETPAR(*)      ,
     &                      SNODE(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             KCOL
C
C***********************************************************************
C
C       --------------------------------------------
C       COMPUTE THE FUNDAMENTAL SUPERNODE PARTITION.
C       --------------------------------------------
        NSUPER = 1
        SNODE(1) = 1
        NOFSUB = COLCNT(1)
        DO  300  KCOL = 2, NEQNS
            IF  ( ETPAR(KCOL-1) .EQ. KCOL )  THEN
                IF  ( COLCNT(KCOL-1) .EQ. COLCNT(KCOL)+1 )  THEN
                    SNODE(KCOL) = NSUPER
                    GO TO 300
                ENDIF
            ENDIF
            NSUPER = NSUPER + 1
            SNODE(KCOL) = NSUPER
            NOFSUB = NOFSUB + COLCNT(KCOL)
  300   CONTINUE
C
      RETURN
      END
