#! /bin/sh -f
#
#   pcxarch.sh - Returns the machine's PCx environmental variable, PCx_ARCH.
#
# This file was copied from the PETSc package
#   http://www.mcs.anl.gov/petsc/petsc.html

if [ -e /bin/uname ]; then
    LARCH=`/bin/uname -s`
    if [ "$LARCH" == "AIX" ]; then
        LARCH="rs6000"
    elif [ "$LARCH" == "HP-UX" ]; then
	LARCH="hpux"
    elif [ "$LARCH" == "Linux" ]; then
	LARCH="linux"
    elif [ "$LARCH" == "FreeBSD" ]; then
	LARCH="freebsd"    
    elif [ "$LARCH" != "IRIX" && "$LARCH" != "IRIX64" ]; then 
        LARCH=`/bin/uname -m`
    fi
elif [ -e /usr/bin/uname ]; then
    LARCH=`/usr/bin/uname`
    if ("$LARCH" == "FreeBSD") then
      LARCH="freebsd"
    else
      echo "Unable to determine machines architecture"
      LARCH="unknown"
  fi
else
    echo "Unable to determine machines architecture"
    LARCH="unknown"
fi

SunOSTest=`expr "$LARCH" : "\(....\)"`
if [ "$SunOSTest" == "sun4" ]; then
  LARCH=sun4
  Version=`/bin/uname -r`
  MajorVersion=`expr "$Version" : "\(.\)"`
  if [ "$MajorVersion" ==  5 ]; then
    LARCH="solaris"
  fi
elif [ "$LARCH" == "AIX" ]; then
   LARCH="rs6000"
elif [ "$LARCH" == "RIOS" ]; then
   LARCH="rs6000"
elif [ "$LARCH" == "sun4m" ]; then
  LARCH="sun4"
elif [ "$LARCH" == "iris4d" ]; then
  LARCH="IRIX"
elif [ "$LARCH" == "Linux" ]; then 
  LARCH="linux"
elif [ "$LARCH" == "CRAY Y-MP" ]; then 
  LARCH="t3d"
fi
#
echo $LARCH
exit 0

