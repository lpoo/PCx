C***********************************************************************
C***********************************************************************
C
C   Written:        October 6, 1996 by SJW. Based on routine BLKSLV of
C                   Esmond G. Ng and Barry W. Peyton.
C
C   Modified:       Sept 30, 1999 to improve efficiency in the case
C                   in which the right-hand side and solution are both
C                   expected to be sparse. Happens a lot in "dense"
C                   column handling.
C
C***********************************************************************
C***********************************************************************
C*********     BLKSLB ... BACK TRIANGULAR SUBSTITUTION        **********
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       GIVEN THE CHOLESKY FACTORIZATION OF A SPARSE SYMMETRIC
C       POSITIVE DEFINITE MATRIX, THIS SUBROUTINE PERFORMS THE
C       BACKWARD TRIANGULAR SUBSTITUTION.  IT USES OUTPUT FROM BLKFCT.
C
C   INPUT PARAMETERS:
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   SUPERNODE PARTITION.
C       (XLINDX,LINDX)  -   ROW INDICES FOR EACH SUPERNODE.
C       (XLNZ,LNZ)      -   CHOLESKY FACTOR.
C
C   UPDATED PARAMETERS:
C       RHS             -   ON INPUT, CONTAINS THE RIGHT HAND SIDE.  ON
C                           OUTPUT, CONTAINS THE SOLUTION.
C
C***********************************************************************
C
      SUBROUTINE  BLKSLB (  NSUPER, XSUPER, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , RHS                             )
C
C***********************************************************************
C
        INTEGER             NSUPER
        INTEGER             LINDX(*)      , XSUPER(*)
        INTEGER             XLINDX(*)     , XLNZ(*)
        DOUBLE PRECISION    LNZ(*)        , RHS(*)
C
C***********************************************************************
C
        INTEGER             FJCOL , I     , IPNT  , IX    , IXSTOP,
     &                      IXSTRT, JCOL  , JPNT  , JSUP  , LJCOL
        DOUBLE PRECISION    T
C
C***********************************************************************
C
        IF  ( NSUPER .LE. 0 )  RETURN
C       -------------------------
C       BACKWARD SUBSTITUTION ...
C       -------------------------
        LJCOL = XSUPER(NSUPER+1) - 1
        DO  600  JSUP = NSUPER, 1, -1
            FJCOL  = XSUPER(JSUP)
            IXSTOP = XLNZ(LJCOL+1) - 1
            JPNT   = XLINDX(JSUP) + (LJCOL - FJCOL)
            DO  500  JCOL = LJCOL, FJCOL, -1
                IXSTRT = XLNZ(JCOL)
                IPNT   = JPNT + 1
                T      = RHS(JCOL)
CDIR$           IVDEP
                DO  400  IX = IXSTRT+1, IXSTOP
                   I = LINDX(IPNT)
                   IF(RHS(I) .NE. 0.D0) T = T - LNZ(IX)*RHS(I)
                   IPNT = IPNT + 1
 400            CONTINUE

                IF(T .NE. 0.D0) THEN
                   RHS(JCOL) = T/LNZ(IXSTRT)
                ELSE
                   RHS(JCOL) = 0.D0
                ENDIF

                IXSTOP    = IXSTRT - 1
                JPNT      = JPNT - 1
 500            CONTINUE

            LJCOL = FJCOL - 1
  600   CONTINUE
C
        RETURN
      END
