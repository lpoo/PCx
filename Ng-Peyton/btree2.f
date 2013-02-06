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
C******     BTREE2 ..... BINARY TREE REPRESENTATION OF ETREE     *******
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       TO DETERMINE A BINARY TREE REPRESENTATION OF THE ELIMINATION 
C       TREE, FOR WHICH EVERY "LAST CHILD" HAS THE MAXIMUM POSSIBLE
C       COLUMN NONZERO COUNT IN THE FACTOR.  THE RETURNED REPRESENTATION 
C       WILL BE GIVEN BY THE FIRST-SON AND BROTHER VECTORS.  THE ROOT OF 
C       THE BINARY TREE IS ALWAYS NEQNS.
C 
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
C                           IT IS ASSUMED THAT PARENT(I) > I EXCEPT OF
C                           THE ROOTS.
C       COLCNT          -   COLUMN NONZERO COUNTS OF THE FACTOR.
C
C   OUTPUT PARAMETERS:
C       FSON            -   THE FIRST SON VECTOR.
C       BROTHR          -   THE BROTHER VECTOR.
C
C   WORKING PARAMETERS:
C       LSON            -   LAST SON VECTOR.
C
C***********************************************************************
C
      SUBROUTINE  BTREE2 (  NEQNS , PARENT, COLCNT, FSON  , BROTHR,
     &                      LSON    )
C
C***********************************************************************
C
        INTEGER             BROTHR(*)     , COLCNT(*)     ,
     &                      FSON(*)       , LSON(*)       ,
     &                      PARENT(*)
C
        INTEGER             NEQNS
C
C***********************************************************************
C
        INTEGER*4           LROOT , NODE  , NDLSON, NDPAR
C
C***********************************************************************
C
        IF  ( NEQNS .LE. 0 )  RETURN
C
        DO  100  NODE = 1, NEQNS
            FSON(NODE) = 0
            BROTHR(NODE) = 0
            LSON(NODE) = 0
  100   CONTINUE
        LROOT = NEQNS
C       ------------------------------------------------------------
C       FOR EACH NODE := NEQNS-1 STEP -1 DOWNTO 1, DO THE FOLLOWING.
C       ------------------------------------------------------------
        IF  ( NEQNS .LE. 1 )  RETURN
        DO  300  NODE = NEQNS-1, 1, -1
            NDPAR = PARENT(NODE)
            IF  ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )  THEN
C               -------------------------------------------------
C               NODE HAS NO PARENT.  GIVEN STRUCTURE IS A FOREST.
C               SET NODE TO BE ONE OF THE ROOTS OF THE TREES.
C               -------------------------------------------------
                BROTHR(LROOT) = NODE
                LROOT = NODE
            ELSE
C               -------------------------------------------
C               OTHERWISE, BECOMES FIRST SON OF ITS PARENT.
C               -------------------------------------------
                NDLSON = LSON(NDPAR)
                IF  ( NDLSON .NE. 0 )  THEN
                    IF  ( COLCNT(NODE) .GE. COLCNT(NDLSON) )  THEN
                        BROTHR(NODE) = FSON(NDPAR)
                        FSON(NDPAR) = NODE
                    ELSE
                        BROTHR(NDLSON) = NODE
                        LSON(NDPAR) = NODE
                    ENDIF
                ELSE
                    FSON(NDPAR) = NODE
                    LSON(NDPAR) = NODE
                ENDIF
            ENDIF
  300   CONTINUE
        BROTHR(LROOT) = 0
C
        RETURN
      END
