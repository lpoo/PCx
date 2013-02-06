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
C******     BFINIT ..... INITIALIZATION FOR BLOCK FACTORIZATION   ******
C***********************************************************************
C***********************************************************************
C
C   PURPOSE:
C       THIS SUBROUTINE COMPUTES ITEMS NEEDED BY THE LEFT-LOOKING
C       BLOCK-TO-BLOCK CHOLESKY FACTORITZATION ROUTINE BLKFCT.
C
C   INPUT PARAMETERS:
C       NEQNS           -   NUMBER OF EQUATIONS.
C       NSUPER          -   NUMBER OF SUPERNODES.
C       XSUPER          -   INTEGER ARRAY OF SIZE (NSUPER+1) CONTAINING
C                           THE SUPERNODE PARTITIONING.
C       SNODE           -   SUPERNODE MEMBERSHIP.
C       (XLINDX,LINDX)  -   ARRAYS DESCRIBING THE SUPERNODAL STRUCTURE.
C       CACHSZ          -   CACHE SIZE (IN KBYTES).
C
C   OUTPUT PARAMETERS:
C       TMPSIZ          -   SIZE OF WORKING STORAGE REQUIRED BY BLKFCT.
C       SPLIT           -   SPLITTING OF SUPERNODES SO THAT THEY FIT
C                           INTO CACHE.
C
C***********************************************************************
C
      SUBROUTINE  BFINIT  ( NEQNS , NSUPER, XSUPER, SNODE , XLINDX, 
     &                      LINDX , CACHSZ, TMPSIZ, SPLIT           )
C
C***********************************************************************
C
        INTEGER     CACHSZ, NEQNS , NSUPER, TMPSIZ
        INTEGER     XLINDX(*)       , XSUPER(*)
        INTEGER     LINDX (*)       , SNODE (*)   ,
     &              SPLIT(*)
C
C***********************************************************************
C
C       ---------------------------------------------------
C       DETERMINE FLOATING POINT WORKING SPACE REQUIREMENT.
C       ---------------------------------------------------
        CALL  FNTSIZ (  NSUPER, XSUPER, SNODE , XLINDX, LINDX ,
     &                  TMPSIZ                                  )
C
C       -------------------------------
C       PARTITION SUPERNODES FOR CACHE.
C       -------------------------------
        CALL  FNSPLT (  NEQNS , NSUPER, XSUPER, XLINDX, CACHSZ,
     &                  SPLIT                                   )
C
        RETURN
      END
