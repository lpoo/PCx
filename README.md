PCx depends of the below programs:

* Fotran compiler, e.g., gfortran >= 4.8.1
* C compiler, e.g., gcc >= 4.8.1
* Make, e.g., GNU Make >= 3.82
* Autoconf, e.g., GNU Autoconf >= 2.69

The PCx's documentation depends of the below programs:

* Sphinx >= 1.2b1
* Doxygen >= 1.8.5
* Breathe >= 1.0.0

To install PCx on a Unix system, follow the below instructions:

    $ git clone https://github.com/r-gaia-cs/PCx.git
    $ ./configure
    $ make

You probably want to take a look at the options of `configure`:

    $ ./configure --help

since some features for debug and profile can be enable using options or
environment variables.

If you got a error like

    configure: error: cannot find install-sh, install.sh, or shtool in "." "./.." "./../.."

try

    $ libtoolize --force
    $ aclocal
    $ autoreconf -f -i -Wall,no-obsolete
    $ ./configure
    $ make

And for the documentation:

    $ make doc

If the building process have finish without error the file `PCx` must be
create. Lets try it with the smallest mps file, `afiro.mps`.

    $ ./PCx afiro

For recompile the source is recommended to remove all the object files.

    $ make clean

To compile the MATLAB interface, follow the below instructions:

    $ make mex

Add PCx/mex to your MATLABPATH environment. Then, at the MATLAB prompt,
type 'help PCx' for syntax information. The syntax is virtually
identical to that of Matlab's linprog (a.k.a. Yin Zhang's LIPSOL).

You can test the interface on the supplied datafiles `afiro.mat`.

=======================================================================

We have tested this procedure on the following architectures:

sun4-class workstations running SunOS;
UltraSparc workstation running Solaris 2.x;
SGI workstations running IRIX 5.3 and 6.4;
HP-9000 workstations;
Pentium PC running Linux;
IBM RS6000 workstations.

Unfortunately, the build procedure is site-specific as well as
architecture-specific, particularly since it involves linking of
Fortran and C source files, so there are no guarantees that it will
work on your machine, even if you have one of the systems mentioned
above. If the "build" fails on your machine, you can download an
executable from the PCx web site for each of the architectures above.
See the PCx home page for details.

Some linkers might issue a warning that there are two main-functions,
chances are that this warning can be safely ignored and the executable
will be fully functional.

If your compiler does not compile the PCx code correctly, try removing
the reference to nullmain.o in the file SRC/Makefile.  Some compilers
seem need this null routine in the linking process.

If you are having difficulty compiling the file timers.c, it is
possible to turn off the timing routines.  Simply change the file
SRC/Makefile.  Change the line

    CFLAGS = -O -D$(PCx_ARCH)

to

    CFLAGS = -O -D$(PCx_ARCH) -DNO_TIMING

To test PCx on one of the input files in the subdirectory "mps/", type

    $ PCx afiro

A User Manual is available in the subdirectory DOC.
