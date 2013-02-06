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
C
      REAL FUNCTION  GTIMER ()
C       --------------------------
C       FOR IBM RS/6000 ...
C       INTEGER     MCLOCK
C       GTIMER = MCLOCK()/100.0
C       --------------------------
C       FOR MOST BERKELEY UNIX ...
        REAL        ETIME
        REAL        VEC(2)
        GTIMER = ETIME(VEC)
C       --------------------------
C       FOR CRAY ...
C       REAL        SECOND
C       GTIMER = SECOND()
C       --------------------------
        RETURN
      END
