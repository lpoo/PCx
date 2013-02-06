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
C****************    FSUP2  ..... FIND SUPERNODES #2   *****************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE IS THE SECOND OF TWO ROUTINES FOR FINDING A
C       MAXIMAL SUPERNODE PARTITION.  IT'S SOLE PURPOSE IS TO 
C       CONSTRUCT THE NEEDED VECTOR OF LENGTH NSUPER: XSUPER(*).  THE
C       FIRST ROUTINE FSUP1 COMPUTES THE NUMBER OF SUPERNODES AND THE 
C       SUPERNODE MEMBERSHIP VECTOR SNODE(*), WHICH IS OF LENGTH NEQNS.
C
C
C   ASSUMPTIONS:
C       THIS ROUTINE ASSUMES A POSTORDERING OF THE ELIMINATION TREE.  IT
C       ALSO ASSUMES THAT THE OUTPUT FROM FSUP1 IS AVAILABLE.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS.
C       (I) NSUPER      -   NUMBER OF SUPERNODES (<= NEQNS).
C       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
C       (I) SNODE(*)    -   ARRAY OF LENGTH NEQNS FOR RECORDING
C                           SUPERNODE MEMBERSHIP.
C
C   OUTPUT PARAMETERS:
C       (I) XSUPER(*)   -   ARRAY OF LENGTH NSUPER+1, CONTAINING THE
C                           SUPERNODE PARTITIONING.
C
C   FIRST CREATED ON    JANUARY 18, 1992.
C   LAST UPDATED ON     NOVEMEBER 22, 1994.
C
C***********************************************************************
C
      SUBROUTINE  FSUP2  (  NEQNS , NSUPER, ETPAR , SNODE , XSUPER  )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             NEQNS , NSUPER
        INTEGER             ETPAR(*)      , SNODE(*)      , 
     &                      XSUPER(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             KCOL  , KSUP  , LSTSUP
C
C***********************************************************************
C
C       -------------------------------------------------
C       COMPUTE THE SUPERNODE PARTITION VECTOR XSUPER(*).
C       -------------------------------------------------
        LSTSUP = NSUPER + 1
        DO  100  KCOL = NEQNS, 1, -1
            KSUP = SNODE(KCOL)
            IF  ( KSUP .NE. LSTSUP )  THEN
                XSUPER(LSTSUP) = KCOL + 1
            ENDIF
            LSTSUP = KSUP
  100   CONTINUE
        XSUPER(1) = 1
C
      RETURN
      END
