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
C**********     CHORDR ..... CHILD REORDERING                ***********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       REARRANGE THE CHILDREN OF EACH VERTEX SO THAT THE LAST ONE 
C       MAXIMIZES (AMONG THE CHILDREN) THE NUMBER OF NONZEROS IN THE 
C       CORRESPONDING COLUMN OF L.  ALSO DETERMINE AN NEW POSTORDERING 
C       BASED ON THE STRUCTURE OF THE MODIFIED ELIMINATION TREE.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
C
C   UPDATED PARAMETERS:
C       (PERM,INVP)     -   ON INPUT, THE GIVEN PERM AND INVERSE PERM
C                           VECTORS.  ON OUTPUT, THE NEW PERM AND
C                           INVERSE PERM VECTORS OF THE NEW
C                           POSTORDERING.
C       COLCNT          -   COLUMN COUNTS IN L UNDER INITIAL ORDERING;
C                           MODIFIED TO REFLECT THE NEW ORDERING.
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
C       BTREE2, EPOST2, INVINV.
C
C***********************************************************************
C
      SUBROUTINE  CHORDR (  NEQNS , XADJ  , ADJNCY, PERM  , INVP  ,
     &                      COLCNT, PARENT, FSON  , BROTHR, INVPOS  )
C
C***********************************************************************
C
        INTEGER             ADJNCY(*)     , BROTHR(*)     ,
     &                      COLCNT(*)     , FSON(*)       , 
     &                      INVP(*)       , INVPOS(*)     , 
     &                      PARENT(*)     , PERM(*)
C
        INTEGER             XADJ(*)
        INTEGER             NEQNS
C
C***********************************************************************
C
C       ----------------------------------------------------------
C       COMPUTE A BINARY REPRESENTATION OF THE ELIMINATION TREE, 
C       SO THAT EACH "LAST CHILD" MAXIMIZES AMONG ITS SIBLINGS THE 
C       NUMBER OF NONZEROS IN THE CORRESPONDING COLUMNS OF L.
C       ----------------------------------------------------------
        CALL  BTREE2  ( NEQNS , PARENT, COLCNT, FSON  , BROTHR, 
     &                  INVPOS                                  )
C
C       ----------------------------------------------------
C       POSTORDER THE ELIMINATION TREE (USING THE NEW BINARY  
C       REPRESENTATION.  
C       ----------------------------------------------------
        CALL  EPOST2  ( NEQNS , FSON  , BROTHR, INVPOS, PARENT, 
     &                  COLCNT, PERM                            ) 
C
C       --------------------------------------------------------
C       COMPOSE THE ORIGINAL ORDERING WITH THE NEW POSTORDERING.
C       --------------------------------------------------------
        CALL  INVINV  ( NEQNS , INVP  , INVPOS, PERM    )
C
        RETURN
      END
