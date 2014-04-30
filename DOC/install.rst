Installing PCx
==============

Executables of PCx for Unix-like environmments (e.g., GNU/Linux) can be built
from source like any other Unix-like program. See the following procedure for
more informations.

#. Clone the git repository. ::

    $ cd $HOME
    $ git clone https://github.com/r-gaia-cs/PCx.git

#. To create the executable ``PCx`` that uses the default NG-Peyton solver ::

    $ cd PCx
    $ ./configure
    $ make

To test PCx with one of the input files in the directory ``mps``, modify the
sample specifications file ``PCx.specs`` if desired and them type ::

    $ ./PCx netlib-afiro

Compilers
---------

The table bellow show the compilers that we already test use to build PCx.

+---------+-----------+---------------+------------------+------------------+------------------+
| Commit  | Processor | C compiler    | Fortran compiler | Compiler success | Solver afiro     |
+---------+-----------+---------------+------------------+------------------+------------------+
| 3fa901d | i386      | gcc-4.4 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.6 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.7 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.8 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         +-----------+---------------+------------------+------------------+------------------+
|         | x86-64    | gcc-4.4 (*)   | gfortran-4.4 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.6 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.7 (*)   | gfortran-4.8 (*) | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.8 (*)   | gfortran-4.8 (*) | YES              | YES              |
+---------+-----------+---------------+------------------+------------------+------------------+
| 5ff4e38 | i386      | clang-3.3 (*) | Not need         | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.4 (*)   | Not need         | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.6 (*)   | Not need         | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.7 (*)   | Not need         | YES              | YES              |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.8 (*)   | Not need         | YES              | YES              |
|         +-----------+---------------+------------------+------------------+------------------+
|         | x86-64    | clang-3.3 (*) | Not need         | YES              | NO               |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.4 (*)   | Not need         | YES              | NO               |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.6 (*)   | Not need         | YES              | NO               |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.7 (*)   | Not need         | YES              | NO               |
|         |           +---------------+------------------+------------------+------------------+
|         |           | gcc-4.8 (*)   | Not need         | YES              | NO               |
+---------+-----------+---------------+------------------+------------------+------------------+

(*) Test in a virtual machine running Debian.
