Installing PCx
==============

Executables for PCx for Unix-like environmments (e.g., GNU/Linux) can be built
from source via the following procedure.

.. note::

   Some environment variables have to be assign::

    $ export FC=gfortran
    $ export CC=gcc

   It will be better if you add the following lines to your ``.bash_profile``::

    export FC=gfortran
    export CC=gcc

#. Clone the git repository. ::

    $ cd $HOME
    $ git clone https://github.com/r-gaia-cs/PCx.git

#. To create the executable ``PCx`` that uses the default NG-Peyton solver ::

    $ cd PCx
    $ ./build.sh

   Because of architectural and environmental differences, it is necessary to
   have a slightly different compiltion procedure for each machine. The
   ``build.sh`` script defines an environment variable ``PCx_ARCH`` and assigns
   it a value to indicate the architecture. ``build.sh`` the invokes the
   ``make`` procedure, with architecture-dependent portions of the makefile
   being retrieved from the subdirectory ``MAKEARCH``. Since the variable
   ``PCx_ARCH`` must be defined for compiling, one should always use ``build``
   instead of ``make`` to compile the program.

To test PCx it on one of the input files int the directory ``mps``, modify the
sample specifications file ``PCx.specs`` if desired, the type ::

    $ ./PCX afiro

or ::

    $ ./PCx 25fc47
