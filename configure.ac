# information on the package
AC_INIT([PCx], [1.1],,[http://pages.cs.wisc.edu/~swright/PCx/])

# Uncompress files in mps folder
echo "Uncompress files in mps folder. It will take some minutes."
cd mps
for f in `ls *.bz2`; do
  if test ! -e ${f/\.bz2/}; then
      echo "    Uncompress $f."
      bunzip2 -c $f > ${f/\.bz2/}
  fi
done
cd ..

# Get the machine's PCx environmental variable, PCx_ARCH.
#
# It was copied from the PETSc package
#   http://www.mcs.anl.gov/petsc/petsc.html
if test -e /bin/uname; then
  PCx_ARCH=`/bin/uname -s`
  if test "$PCx_ARCH" == "AIX"; then
    PCx_ARCH="rs6000"
  elif test "$PCx_ARCH" == "HP-UX"; then
    PCx_ARCH="hpux"
  elif test "$PCx_ARCH" == "Linux"; then
    PCx_ARCH="linux"
  elif test "$PCx_ARCH" == "FreeBSD"; then
    PCx_ARCH="freebsd"
  elif test "$PCx_ARCH" != "IRIX" && "$PCx_ARCH" != "IRIX64"; then
    PCx_ARCH=`/bin/uname -m`
  fi
elif test -e /usr/bin/uname; then
  PCx_ARCH=`/usr/bin/uname`
  if test "$PCx_ARCH" == "FreeBSD"; then
    PCx_ARCH="freebsd"
  else
    echo "Unable to determine machines architecture"
    PCx_ARCH="unknown"
  fi
else
  echo "Unable to determine machines architecture"
  PCx_ARCH="unknown"
fi

SunOSTest=`expr "$PCx_ARCH" : "\(....\)"`
if test "$SunOSTest" == "sun4"; then
  PCx_ARCH=sun4
  Version=`/bin/uname -r`
  MajorVersion=`expr "$Version" : "\(.\)"`
  if test "$MajorVersion" ==  5; then
    PCx_ARCH="solaris"
  fi
elif test "$PCx_ARCH" == "AIX"; then
   PCx_ARCH="rs6000"
elif test "$PCx_ARCH" == "RIOS"; then
   PCx_ARCH="rs6000"
elif test "$PCx_ARCH" == "sun4m"; then
  PCx_ARCH="sun4"
elif test "$PCx_ARCH" == "iris4d"; then
  PCx_ARCH="IRIX"
elif test "$PCx_ARCH" == "Linux"; then
  PCx_ARCH="linux"
elif test "$PCx_ARCH" == "CRAY Y-MP"; then
  PCx_ARCH="t3d"
fi
AC_SUBST(PCx_ARCH)

# Configure the Cholesky solver
echo 'Cheking for IBM WSSMP Cholesky solver';
if test -z "$WSSMP_LIB"; then
  if test -f wssmp/libwssmpp2.a; then
    WSSMP_LIB=../wssmp/libwssmpp2.a
    echo 'WSSMP library found. Using `../wssmp/libwssmpp2.a`';
  else
    if test -f wssmp/libwssmp.a; then
      WSSMP_LIB=../wssmp/libwssmp.a;
      echo 'WSSMP library found. Using `../wssmp/libwssmp.a`';
    else
      echo 'WSSMP library not found in ./wssmp/';
    fi
  fi
fi
AC_SUBST(WSSMP_LIB)

echo 'Checking for Ng Peyton Cholesky solver';
if test -z "$NG_LIB"; then
  if test -d Ng-Peyton; then
    NG_LIB=../Ng-Peyton/cholesky.a
    echo "Ng-Peyton found";
  else
    echo "No directory ./Ng-Peyton";
  fi
fi
AC_SUBST(NG_LIB)

echo 'Checking for user-supplied solver';
if test ! -f ./SRC/mysolver.c; then
  echo 'File ./SRC/mysolver.c not present';
fi
if test ! -f ./mysolver/libmysolver.a; then
  echo 'Library ./mysolver/libmysolver.a not present';
fi

echo 'Checking for the Matlab Interface';
if test -z "$NG_LIB" ; then
  echo "The Matlab Interface depends of Ng Peyton Cholesky solver";
fi


# Configure target dir
TARGETDIR=.
AC_SUBST(TARGETDIR)

# Library
AC_CONFIG_MACRO_DIR([m4])
AM_DISABLE_SHARED
AM_PROG_LIBTOOL

# More options
AC_ARG_ENABLE([warnings],
    AS_HELP_STRING([--disable-warnings], [Disable warnings]),
    [CFLAGS="";
     FCFLAGS=""],
    [CFLAGS="-Wall ";
     FCFLAGS="-Wall "])
AC_ARG_ENABLE([debug],
    AS_HELP_STRING([--enable-debug], [Enable debug features]),
    [CFLAGS+="-g -O0 -DDEBUG ";
     FCFLAGS+="-g -O0 "],
    [CFLAGS+="-O3 ";
     FCFLAGS+="-O3 "])
AC_ARG_ENABLE([profile],
    AS_HELP_STRING([--enable-profile], [Enable profile features]),
    [CFLAGS+="-pg ";
     FCFLAGS+="-pg "],
    [])
AC_ARG_ENABLE([timing],
    AS_HELP_STRING([--enable-timing], [Enable full native timing profile in log file]),
    [CFLAGS+="-DTIMING_PROFILE "],
    [])
AC_ARG_ENABLE([verbose],
    AS_HELP_STRING([--enable-verbose], [Enable verbose features]),
    [CFLAGS+="-DPREPROCESS_VERBOSE "],
    [])
AC_ARG_VAR(DECIMAL_WIDTH,
   Width of the decimal part of time in log (no more than 6))
DECIMAL_WIDTH=${DECIMAL_WIDTH-6}
echo $DECIMAL_WIDTH
CFLAGS+="-DDECIMAL_WIDTH=$DECIMAL_WIDTH "
#AC_SUBST(DECIMAL_WIDTH)

# Checks for programs
echo "Testing for a C compiler"
# Since we need '-z muldefs' we will prefer gcc
AC_PROG_CC([gcc])
echo "Testing for a FORTRAN compiler"
# We need gfortran since f77 didn't support '-z muldefs'
AC_PROG_FC([gfortran])

# Checks for libraries
# Checks for header files
# Checks for types
# Checks for structures
# Checks for compiler characteristics
# Checks for library functions
# Checks for system
# Services

AC_CONFIG_FILES([Makefile Ng-Peyton/Makefile SRC/Makefile])

AC_OUTPUT
