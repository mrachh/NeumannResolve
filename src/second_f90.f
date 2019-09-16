      REAL *8             FUNCTION SECOND( )
*
*  Purpose
*  =======
*
*  SECOND returns the wall time for a process in seconds.
*  This version gets the time from the system function SYSTEM_CLOCK.
*
* =====================================================================
*
        INTEGER ICOUNT, ICOUNT_RATE, ICOUNT_MAX
c     
        call system_clock(icount,icount_rate,icount_max)
        second = icount/dble(icount_rate)
c
c
      RETURN
*
*     End of SECOND
*
      END
