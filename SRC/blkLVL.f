C Interface to the BLKFCT routine.

C PCx beta-0.1   5/1/96
C 
C Authors: Joe Czyzyk, Sanjay Mehrotra, Steve Wright.
C
C (C) 1996 University of Chicago. See COPYRIGHT in main directory.

C It's called from the C routine with all the required
C parameters except for the matrix-matrix and matrix-vector
C multiply routines MMPYN and SMXPY. Instead, there's a 
C flag "LEVEL" which indicates what level of loop unrolling
C is required. This routine plugs in the appropriate two
C routines and then calls BLKFCT

      SUBROUTINE  BLKLVL (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT , 
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , LEVEL   )
C
C***********************************************************************
C
C       -----------
C       PARAMETERS.
C       -----------
C
        INTEGER             XLINDX(*)     , XLNZ(*)
        INTEGER             IWORK(*)      , LINDX(*)      , 
     &                      SNODE(*)      , SPLIT(*)      , 
     &                      XSUPER(*)
        INTEGER             IFLAG , IWSIZ , NEQNS , NSUPER, TMPSIZ
        INTEGER             LEVEL
        DOUBLE PRECISION    LNZ(*)        , TMPVEC(*)

C***********************************************************************

        EXTERNAL            SMXPY1, SMXPY2, SMXPY4, SMXPY8
        EXTERNAL            MMPY1 , MMPY2 , MMPY4 , MMPY8

        IF  ( LEVEL .NE. 1  .AND.  LEVEL .NE. 2  .AND.
     &        LEVEL .NE. 4  .AND.  LEVEL .NE. 8        )  THEN
            WRITE (*,23)  LEVEL
   23       FORMAT ('*** LOOP UNROLLING LEVEL = ', I5,' ***'/,
     &              '*** SHOULD HAVE LEVEL = 1, 2, 4, OR 8 ***' )
            STOP
        ENDIF

        IF  ( LEVEL .EQ. 1 )  THEN
            CALL  BLKFCT (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT ,
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY1 ,
     &                      SMXPY1                                  )
        ELSEIF  ( LEVEL .EQ. 2 )  THEN
            CALL  BLKFCT (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT ,
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY2 ,
     &                      SMXPY2                                  )
        ELSEIF  ( LEVEL .EQ. 4 )  THEN
            CALL  BLKFCT (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT ,
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY4 ,
     &                      SMXPY4                                  )
        ELSEIF  ( LEVEL .EQ. 8 )  THEN
            CALL  BLKFCT (  NEQNS , NSUPER, XSUPER, SNODE , SPLIT ,
     &                      XLINDX, LINDX , XLNZ  , LNZ   , IWSIZ ,
     &                      IWORK , TMPSIZ, TMPVEC, IFLAG , MMPY8 ,
     &                      SMXPY8                                  )
        ENDIF

        RETURN
        END
