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
C***************     ETPOST ..... ETREE POSTORDERING     ***************
C***********************************************************************
C***********************************************************************
C
C   WRITTEN BY JOSEPH LIU (SEPT 17, 1986)
C
C   PURPOSE:
C       BASED ON THE BINARY REPRESENTATION (FIRST-SON,BROTHER) OF
C       THE ELIMINATION TREE, A POSTORDERING IS DETERMINED. THE
C       CORRESPONDING PARENT VECTOR IS ALSO MODIFIED TO REFLECT
C       THE REORDERING.
C
C   INPUT PARAMETERS:
C       ROOT            -   ROOT OF THE ELIMINATION TREE (USUALLY IT
C                           IS NEQNS).
C       FSON            -   THE FIRST SON VECTOR.
C       BROTHR          -   THE BROTHR VECTOR.
C
C   UPDATED PARAMETERS:
C       PARENT          -   THE PARENT VECTOR.
C
C   OUTPUT PARAMETERS:
C       INVPOS          -   INVERSE PERMUTATION FOR THE POSTORDERING.
C
C   WORKING PARAMETERS:
C       STACK           -   THE STACK FOR POSTORDER TRAVERSAL OF THE
C                           TREE.
C
C***********************************************************************
C
      SUBROUTINE  ETPOST (  ROOT  , FSON  , BROTHR, INVPOS, PARENT,
     &                      STACK                                   )
C
C***********************************************************************
C
        INTEGER*4           BROTHR(*)     , FSON(*)       ,
     &                      INVPOS(*)     , PARENT(*)     ,
     &                      STACK(*)
C
        INTEGER*4           ROOT
C
C***********************************************************************
C
        INTEGER*4           ITOP  , NDPAR , NODE  , NUM   , NUNODE
C
C***********************************************************************
C
        NUM = 0
        ITOP = 0
        NODE = ROOT
C       -------------------------------------------------------------
C       TRAVERSE ALONG THE FIRST SONS POINTER AND PUSH THE TREE NODES
C       ALONG THE TRAVERSAL INTO THE STACK.
C       -------------------------------------------------------------
  100   CONTINUE
            ITOP = ITOP + 1
            STACK(ITOP) = NODE
            NODE = FSON(NODE)
            IF  ( NODE .GT. 0 )  GO TO 100
C           ----------------------------------------------------------
C           IF POSSIBLE, POP A TREE NODE FROM THE STACK AND NUMBER IT.
C           ----------------------------------------------------------
  200       CONTINUE
                IF  ( ITOP .LE. 0 )  GO TO 300
                NODE = STACK(ITOP)
                ITOP = ITOP - 1
                NUM = NUM + 1
                INVPOS(NODE) = NUM
C               ----------------------------------------------------
C               THEN, TRAVERSE TO ITS YOUNGER BROTHER IF IT HAS ONE.
C               ----------------------------------------------------
                NODE = BROTHR(NODE)
                IF  ( NODE .LE. 0 )  GO TO 200
            GO TO 100
C
  300   CONTINUE
C       ------------------------------------------------------------
C       DETERMINE THE NEW PARENT VECTOR OF THE POSTORDERING.  BROTHR
C       IS USED TEMPORARILY FOR THE NEW PARENT VECTOR.
C       ------------------------------------------------------------
        DO  400  NODE = 1, NUM
            NUNODE = INVPOS(NODE)
            NDPAR = PARENT(NODE)
            IF  ( NDPAR .GT. 0 )  NDPAR = INVPOS(NDPAR)
            BROTHR(NUNODE) = NDPAR
  400   CONTINUE
C
        DO  500  NUNODE = 1, NUM
            PARENT(NUNODE) = BROTHR(NUNODE)
  500   CONTINUE
C
        RETURN
      END
