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
C******     FNTSIZ ..... COMPUTE WORK STORAGE SIZE FOR BLKFCT     ******
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE DETERMINES THE SIZE OF THE WORKING STORAGE
C       REQUIRED BY BLKFCT.
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
C                           THE SUPERNODE PARTITIONING.
C       SNODE           -   SUPERNODE MEMBERSHIP.
C       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE.
C
C   OUTPUT PARAMETERS:
C       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT.
C
C***********************************************************************
C
        SUBROUTINE  FNTSIZ ( NSUPER, XSUPER, SNODE , XLINDX, 
     &                       LINDX , TMPSIZ  )
C
C***********************************************************************
C
        INTEGER     NSUPER, TMPSIZ
        INTEGER     XLINDX(*)       , XSUPER(*)
        INTEGER     LINDX (*)       , SNODE (*)
C
        INTEGER     BOUND , CLEN  , CURSUP, I     , IBEGIN, IEND  , 
     &              KSUP  , LENGTH, NCOLS , NXTSUP, 
     &              TSIZE , WIDTH
C
C***********************************************************************
C
C       RETURNS SIZE OF TEMP ARRAY USED BY BLKFCT FACTORIZATION ROUTINE.
C       NOTE THAT THE VALUE RETURNED IS AN ESTIMATE, THOUGH IT IS USUALLY
C       TIGHT.
C
C       ----------------------------------------
C       COMPUTE SIZE OF TEMPORARY STORAGE VECTOR
C       NEEDED BY BLKFCT.
C       ----------------------------------------
        TMPSIZ = 0
        DO  500  KSUP = NSUPER, 1, -1
            NCOLS = XSUPER(KSUP+1) - XSUPER(KSUP)
            IBEGIN = XLINDX(KSUP) + NCOLS
            IEND = XLINDX(KSUP+1) - 1
            LENGTH = IEND - IBEGIN + 1
            BOUND = LENGTH * (LENGTH + 1) / 2
            IF  ( BOUND .GT. TMPSIZ )  THEN
                CURSUP = SNODE(LINDX(IBEGIN))
                CLEN = XLINDX(CURSUP+1) - XLINDX(CURSUP)
                WIDTH = 0
                DO  400  I = IBEGIN, IEND
                    NXTSUP = SNODE(LINDX(I))
                    IF  ( NXTSUP .EQ. CURSUP )  THEN
                        WIDTH = WIDTH + 1
                        IF  ( I .EQ. IEND )  THEN
                            IF  ( CLEN .GT. LENGTH )  THEN
                                TSIZE = LENGTH * WIDTH - 
     &                                  (WIDTH - 1) * WIDTH / 2
                                TMPSIZ = MAX ( TSIZE , TMPSIZ )
                            ENDIF
                        ENDIF
                    ELSE
                        IF  ( CLEN .GT. LENGTH )  THEN
                            TSIZE = LENGTH * WIDTH - 
     &                              (WIDTH - 1) * WIDTH / 2
                            TMPSIZ = MAX ( TSIZE , TMPSIZ )
                        ENDIF
                        LENGTH = LENGTH - WIDTH
                        BOUND = LENGTH * (LENGTH + 1) / 2
                        IF  ( BOUND .LE. TMPSIZ )  GO TO 500
                        WIDTH = 1
                        CURSUP = NXTSUP
                        CLEN = XLINDX(CURSUP+1) - XLINDX(CURSUP)
                    ENDIF
  400           CONTINUE
            ENDIF
  500   CONTINUE
C
        RETURN
        END
