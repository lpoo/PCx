C***********************************************************************
C***********************************************************************
C
C   Version:        0.4
C   Last modified:  December 27, 1994
C   Authors:        Joseph W.H. Liu
C
C   Mathematical Sciences Section, Oak Ridge National Laboratory
C
C***********************************************************************
C***********************************************************************
C**********     ETORDR ..... ELIMINATION TREE REORDERING     ***********
C***********************************************************************
C***********************************************************************
C
C   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
C
C   PURPOSE:
C       TO DETERMINE AN EQUIVALENT REORDERING BASED ON THE STRUCTURE OF
C       THE ELIMINATION TREE.  A POSTORDERING OF THE GIVEN ELIMINATION
C       TREE IS RETURNED.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
C
C   UPDATED PARAMETERS:
C       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM
C                           VECTORS.  ON OUTPUT, THE NEW PERM AND
C                           INVERSE PERM VECTORS OF THE EQUIVALENT
C                           ORDERING.
C
C   OUTPUT PARAMETERS:
C       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE
C                           ASSOCIATED WITH THE NEW ORDERING.
C
C   WORKING PARAMETERS:
C       FSON            -   THE FIRST SON VECTOR.
C       BROTHR          -   THE BROTHER VECTOR.
C       INVPOS          -   THE INVERSE PERM VECTOR FOR THE
C                           POSTORDERING.
C
C   PROGRAM SUBROUTINES:
C       BETREE, ETPOST, ETREE , INVINV.
C
C***********************************************************************
C
      SUBROUTINE  ETORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  ,
     &                      PARENT, FSON  , BROTHR, INVPOS          )
C
C***********************************************************************
C
        INTEGER*4           ADJNCY(*)     , BROTHR(*)     ,
     &                      FSON(*)       , INVP(*)       ,
     &                      INVPOS(*)     , PARENT(*)     ,
     &                      PERM(*)
C
        INTEGER*4           XADJ(*)
        INTEGER*4           NEQNS
C
C***********************************************************************
C
C       -----------------------------
C       COMPUTE THE ELIMINATION TREE.
C       -----------------------------
        CALL  ETREE ( NEQNS, XADJ, ADJNCY, PERM, INVP, PARENT, INVPOS )
C
C       --------------------------------------------------------
C       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE.
C       --------------------------------------------------------
        CALL  BETREE ( NEQNS, PARENT, FSON, BROTHR )
C
C       -------------------------------
C       POSTORDER THE ELIMINATION TREE.
C       -------------------------------
        CALL  ETPOST ( NEQNS, FSON, BROTHR, INVPOS, PARENT, PERM )
C
C       --------------------------------------------------------
C       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING.
C       --------------------------------------------------------
        CALL  INVINV ( NEQNS, INVP, INVPOS, PERM )
C
        RETURN
      END
