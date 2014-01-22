/* timing routines
 *
 * PCx 1.1 11/97
 *
 * Authors: Joe Czyzyk, Sanjay Mehrotra, Michael Wagner, Steve Wright.
 * 
 * (C) 1996 University of Chicago. See COPYRIGHT in main directory.
 */

/*
#define NO_TIMING
*/


#ifdef NO_TIMING

GetTime(User, System)
  double         *User, *System;
{
  *User   = 0.0;
  *System = 0.0;
}

#else             /* perform timing */

#ifndef solaris   /* non-Solaris machines */

#include <stdio.h>
#include <sys/time.h>     /* for getrusage */
#include <sys/resource.h>
#if defined(hpux) || defined(__hpux)
#include <sys/syscall.h>
#define getrusage(a, b)     syscall(SYS_GETRUSAGE, a, b)
#endif

int GetTime(User, System)
  double         *User, *System;
{
  struct rusage   Time;

if (getrusage(RUSAGE_SELF, &Time))
    printf("Error getting time.\n"); 
  *User = Time.ru_utime.tv_sec + Time.ru_utime.tv_usec / 1E6;
  *System = Time.ru_stime.tv_sec + Time.ru_stime.tv_usec / 1E6;
  return 0;
}

#else

/* timer routine for solaris */

#include <stdio.h>

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/signal.h>
#include <sys/fault.h>
#include <sys/syscall.h>

#include <sys/stat.h>
#include <fcntl.h>

/*#ifndef RUSAGE_SELF  dm*/
#include <sys/procfs.h>
/*#endif dm*/

GetTime(User, System)
double *User, *System;
{

  int             fd;
  char            proc[50];
  prusage_t       prusage;
  
  sprintf(proc,"/proc/%d", getpid());

  if ((fd = open(proc,O_RDONLY)) == -1) 
  {
    printf("Error opening process file.\n");
    *User = 0.0;
    *System = 0.0;
    close(fd);
    return;
  }

  if (ioctl(fd, PIOCUSAGE, &prusage) == -1) {
    printf("Error performing ioctl.");
    *User = 0.0;
    *System = 0.0;
    return;
  }

  *User   = prusage.pr_utime.tv_sec + prusage.pr_utime.tv_nsec / 1.0e9;
  *System = prusage.pr_stime.tv_sec + prusage.pr_stime.tv_nsec / 1.0e9;
  close(fd);
}

#endif /* solaris */


#endif /* NO_TIMING */
