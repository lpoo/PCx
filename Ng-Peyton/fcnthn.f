C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  January 12, 1995
C   Authors:        Esmond G. Ng and Barry W. Peyton
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C**************     FCNTHN  ..... FIND NONZERO COUNTS    ***************
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE DETERMINES THE ROW COUNTS AND COLUMN COUNTS IN
C       THE CHOLESKY FACTOR.  IT USES A DISJOINT SET UNION ALGORITHM.
C
C       TECHNIQUES:
C       1) SUPERNODE DETECTION.
C       2) PATH HALVING.
C       3) NO UNION BY RANK.
C
C   NOTES:
C       1) ASSUMES A POSTORDERING OF THE ELIMINATION TREE.
C
C   INPUT PARAMETERS:
C       (I) NEQNS       -   NUMBER OF EQUATIONS.
C       (I) ADJLEN      -   LENGTH OF ADJACENCY STRUCTURE.
C       (I) XADJ(*)     -   ARRAY OF LENGTH NEQNS+1, CONTAINING POINTERS
C                           TO THE ADJACENCY STRUCTURE.
C       (I) ADJNCY(*)   -   ARRAY OF LENGTH XADJ(NEQNS+1)-1, CONTAINING
C                           THE ADJACENCY STRUCTURE.
C       (I) PERM(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           POSTORDERING.
C       (I) INVP(*)     -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           INVERSE OF THE POSTORDERING.
C       (I) ETPAR(*)    -   ARRAY OF LENGTH NEQNS, CONTAINING THE
C                           ELIMINATION TREE OF THE POSTORDERED MATRIX.
C
C   OUTPUT PARAMETERS:
C       (I) ROWCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH ROW OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) COLCNT(*)   -   ARRAY OF LENGTH NEQNS, CONTAINING THE NUMBER
C                           OF NONZEROS IN EACH COLUMN OF THE FACTOR,
C                           INCLUDING THE DIAGONAL ENTRY.
C       (I) NLNZ        -   NUMBER OF NONZEROS IN THE FACTOR, INCLUDING
C                           THE DIAGONAL ENTRIES.
C
C   WORK PARAMETERS:
C       (I) SET(*)      -   ARRAY OF LENGTH NEQNS USED TO MAINTAIN THE
C                           DISJOINT SETS (I.E., SUBTREES).
C       (I) PRVLF(*)    -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
C                           PREVIOUS LEAF OF EACH ROW SUBTREE.
C       (I) LEVEL(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE LEVEL
C                           (DISTANCE FROM THE ROOT).
C       (I) WEIGHT(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING WEIGHTS
C                           USED TO COMPUTE COLUMN COUNTS.
C       (I) FDESC(*)    -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
C                           FIRST (I.E., LOWEST-NUMBERED) DESCENDANT.
C       (I) NCHILD(*)   -   ARRAY OF LENGTH NEQNS+1 CONTAINING THE
C                           NUMBER OF CHILDREN.
C       (I) PRVNBR(*)   -   ARRAY OF LENGTH NEQNS USED TO RECORD THE
C                           PREVIOUS ``LOWER NEIGHBOR'' OF EACH NODE.
C
C   FIRST CREATED ON    APRIL 12, 1990.
C   LAST UPDATED ON     JANUARY 12, 1995.
C
C***********************************************************************
C
      SUBROUTINE FCNTHN  (  NEQNS , ADJLEN, XADJ  , ADJNCY, PERM  ,
     &                      INVP  , ETPAR , ROWCNT, COLCNT, NLNZ  ,
     &                      SET   , PRVLF , LEVEL , WEIGHT, FDESC ,
     &                      NCHILD, PRVNBR                          )
C
C       -----------
C       PARAMETERS.
C       -----------
        INTEGER             ADJLEN, NEQNS , NLNZ
        INTEGER             ADJNCY(ADJLEN)  , COLCNT(NEQNS) ,
     &                      ETPAR(NEQNS)    , FDESC(0:NEQNS),
     &                      INVP(NEQNS)     , LEVEL(0:NEQNS),
     &                      NCHILD(0:NEQNS) , PERM(NEQNS)   ,
     &                      PRVLF(NEQNS)    , PRVNBR(NEQNS) ,
     &                      ROWCNT(NEQNS)   , SET(NEQNS)    ,
     &                      WEIGHT(0:NEQNS)
        INTEGER             XADJ(*)
C
C       ----------------
C       LOCAL VARIABLES.
C       ----------------
        INTEGER             HINBR , IFDESC, J     , JSTOP , JSTRT , 
     &                      K     , LAST1 , LAST2 , LCA   , LFLAG , 
     &                      LOWNBR, OLDNBR, PARENT, PLEAF , TEMP  , 
     &                      XSUP
C
C***********************************************************************
C
C       --------------------------------------------------
C       COMPUTE LEVEL(*), FDESC(*), NCHILD(*).
C       INITIALIZE XSUP, ROWCNT(*), COLCNT(*),
C                  SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*).
C       --------------------------------------------------
        XSUP = 1
        LEVEL(0) = 0
        DO  100  K = NEQNS, 1, -1
            ROWCNT(K) = 1
            COLCNT(K) = 0
            SET(K) = K
            PRVLF(K) = 0
            LEVEL(K) = LEVEL(ETPAR(K)) + 1
            WEIGHT(K) = 1
            FDESC(K) = K
            NCHILD(K) = 0
            PRVNBR(K) = 0
  100   CONTINUE
        NCHILD(0) = 0
        FDESC(0) = 0
        DO  200  K = 1, NEQNS
            PARENT = ETPAR(K)
            WEIGHT(PARENT) = 0
            NCHILD(PARENT) = NCHILD(PARENT) + 1
            IFDESC = FDESC(K)
            IF  ( IFDESC .LT. FDESC(PARENT) )  THEN
                FDESC(PARENT) = IFDESC
            ENDIF
  200   CONTINUE
C       ------------------------------------
C       FOR EACH ``LOW NEIGHBOR'' LOWNBR ...
C       ------------------------------------
        DO  600  LOWNBR = 1, NEQNS
            LFLAG = 0
            IFDESC = FDESC(LOWNBR)
            OLDNBR = PERM(LOWNBR)
            JSTRT = XADJ(OLDNBR)
            JSTOP = XADJ(OLDNBR+1) - 1
C           -----------------------------------------------
C           FOR EACH ``HIGH NEIGHBOR'', HINBR OF LOWNBR ...
C           -----------------------------------------------
            DO  500  J = JSTRT, JSTOP
                HINBR = INVP(ADJNCY(J))
                IF  ( HINBR .GT. LOWNBR )  THEN
                    IF  ( IFDESC .GT. PRVNBR(HINBR) )  THEN
C                       -------------------------
C                       INCREMENT WEIGHT(LOWNBR).
C                       -------------------------
                        WEIGHT(LOWNBR) = WEIGHT(LOWNBR) + 1
                        PLEAF = PRVLF(HINBR)
C                       -----------------------------------------
C                       IF HINBR HAS NO PREVIOUS ``LOW NEIGHBOR'' 
C                       THEN ...
C                       -----------------------------------------
                        IF  ( PLEAF .EQ. 0 )  THEN
C                           -----------------------------------------
C                           ... ACCUMULATE LOWNBR-->HINBR PATH LENGTH 
C                               IN ROWCNT(HINBR).
C                           -----------------------------------------
                            ROWCNT(HINBR) = ROWCNT(HINBR) +
     &                                  LEVEL(LOWNBR) - LEVEL(HINBR)
                        ELSE
C                           -----------------------------------------
C                           ... OTHERWISE, LCA <-- FIND(PLEAF), WHICH 
C                               IS THE LEAST COMMON ANCESTOR OF PLEAF 
C                               AND LOWNBR.
C                               (PATH HALVING.)
C                           -----------------------------------------
                            LAST1 = PLEAF
                            LAST2 = SET(LAST1)
                            LCA = SET(LAST2)
  300                       CONTINUE
                                IF  ( LCA .NE. LAST2 )  THEN
                                    SET(LAST1) = LCA
                                    LAST1 = LCA
                                    LAST2 = SET(LAST1)
                                    LCA = SET(LAST2)
                                    GO TO 300
                                ENDIF
C                           -------------------------------------
C                           ACCUMULATE PLEAF-->LCA PATH LENGTH IN 
C                           ROWCNT(HINBR).
C                           DECREMENT WEIGHT(LCA).
C                           -------------------------------------
                            ROWCNT(HINBR) = ROWCNT(HINBR)
     &                                  + LEVEL(LOWNBR) - LEVEL(LCA)
                            WEIGHT(LCA) = WEIGHT(LCA) - 1
                        ENDIF
C                       ----------------------------------------------
C                       LOWNBR NOW BECOMES ``PREVIOUS LEAF'' OF HINBR.
C                       ----------------------------------------------
                        PRVLF(HINBR) = LOWNBR
                        LFLAG = 1
                    ENDIF
C                   --------------------------------------------------
C                   LOWNBR NOW BECOMES ``PREVIOUS NEIGHBOR'' OF HINBR.
C                   --------------------------------------------------
                    PRVNBR(HINBR) = LOWNBR
                ENDIF
  500       CONTINUE
C           ----------------------------------------------------
C           DECREMENT WEIGHT ( PARENT(LOWNBR) ).
C           SET ( P(LOWNBR) ) <-- SET ( P(LOWNBR) ) + SET(XSUP).
C           ----------------------------------------------------
            PARENT = ETPAR(LOWNBR)
            WEIGHT(PARENT) = WEIGHT(PARENT) - 1
            IF  ( LFLAG .EQ. 1     .OR.
     &            NCHILD(LOWNBR) .GE. 2 )  THEN
                XSUP = LOWNBR
            ENDIF
            SET(XSUP) = PARENT
  600   CONTINUE
C       ---------------------------------------------------------
C       USE WEIGHTS TO COMPUTE COLUMN (AND TOTAL) NONZERO COUNTS.
C       ---------------------------------------------------------
        NLNZ = 0
        DO  700  K = 1, NEQNS
            TEMP = COLCNT(K) + WEIGHT(K)
            COLCNT(K) = TEMP
            NLNZ = NLNZ + TEMP
            PARENT = ETPAR(K)
            IF  ( PARENT .NE. 0 )  THEN
                COLCNT(PARENT) = COLCNT(PARENT) + TEMP
            ENDIF
  700   CONTINUE
C
        RETURN
      END
