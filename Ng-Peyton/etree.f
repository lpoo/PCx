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
C****************     ETREE ..... ELIMINATION TREE     *****************
C***********************************************************************
C***********************************************************************
C
C   WRITTEN BY JOSEPH LIU (JUL 17, 1985)
C
C   PURPOSE:
C       TO DETERMINE THE ELIMINATION TREE FROM A GIVEN ORDERING AND
C       THE ADJACENCY STRUCTURE.  THE PARENT VECTOR IS RETURNED.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
C       (PERM,INVP)     -   PERMUTATION AND INVERSE PERMUTATION VECTORS
C
C   OUTPUT PARAMETERS:
C       PARENT          -   THE PARENT VECTOR OF THE ELIMINATION TREE.
C
C   WORKING PARAMETERS:
C       ANCSTR          -   THE ANCESTOR VECTOR.
C
C***********************************************************************
C
      SUBROUTINE  ETREE (   NEQNS , XADJ  , ADJNCY, PERM  , INVP  ,
     &                      PARENT, ANCSTR                          )
C
C***********************************************************************
C
        INTEGER*4           ADJNCY(*)     , ANCSTR(*)     ,
     &                      INVP(*)       , PARENT(*)     ,
     &                      PERM(*)
C
        INTEGER*4           NEQNS
        INTEGER*4           XADJ(*)
C
C***********************************************************************
C
        INTEGER*4           I     , J     , JSTOP , JSTRT , NBR   ,
     &                      NEXT  , NODE
C
C***********************************************************************
C
        IF  ( NEQNS .LE. 0 )  RETURN
C
        DO  400  I = 1, NEQNS
            PARENT(I) = 0
            ANCSTR(I) = 0
            NODE = PERM(I)
C
            JSTRT = XADJ(NODE)
            JSTOP = XADJ(NODE+1) - 1
            IF  ( JSTRT .LE. JSTOP )  THEN
                DO  300  J = JSTRT, JSTOP
                    NBR = ADJNCY(J)
                    NBR = INVP(NBR)
                    IF  ( NBR .LT. I )  THEN
C                       -------------------------------------------
C                       FOR EACH NBR, FIND THE ROOT OF ITS CURRENT
C                       ELIMINATION TREE.  PERFORM PATH COMPRESSION
C                       AS THE SUBTREE IS TRAVERSED.
C                       -------------------------------------------
  100                   CONTINUE
                            IF  ( ANCSTR(NBR) .EQ. I )  GO TO 300
                            IF  ( ANCSTR(NBR) .GT. 0 )  THEN
                                NEXT = ANCSTR(NBR)
                                ANCSTR(NBR) = I
                                NBR = NEXT
                                GO TO 100
                            ENDIF
C                       --------------------------------------------
C                       NOW, NBR IS THE ROOT OF THE SUBTREE.  MAKE I
C                       THE PARENT NODE OF THIS ROOT.
C                       --------------------------------------------
                        PARENT(NBR) = I
                        ANCSTR(NBR) = I
                    ENDIF
  300           CONTINUE
            ENDIF
  400   CONTINUE
C
        RETURN
      END
