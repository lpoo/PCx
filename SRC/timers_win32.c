// Replaces timers.c - source in timers_win32.c

#define _CRT_SECURE_NO_WARNINGS 1 // Fernando
#define _CRT_NONSTDC_NO_WARNINGS 1 // Fernando

#ifdef WIN32

#include <stdio.h>
#include <sys/timeb.h>
#include <time.h>

int GetTime(User, System)
  double         *User, *System;
{
   struct __timeb64 timebuffer;
    _ftime64( &timebuffer );
	*User = timebuffer.time +  timebuffer.millitm/1000.0;
	*System = 0;
	return 0;
}
#else
GetTime(User, System)
  double         *User, *System;
{
  *User   = 0.0;
  *System = 0.0;
}

#endif
