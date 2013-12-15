Installing PCx
==============

Dependencies
------------

PCx depends of the below programs:

* Fotran compiler, e.g., gfortran
* C compiler, e.g., gcc
* Make, e.g., GNU Make
* Autoconf, e.g., GNU Autoconf

The PCx's documentation depends of the below programs:

* Sphinx
* Doxygen
* Breathe

Basic
-----

To install PCx on a GNU/Linux system, follow the below instructions:

    $ git clone https://github.com/r-gaia-cs/PCx.git
    $ ./configure
    $ make

If you got a error like

    configure: error: cannot find install-sh, install.sh, or shtool in "." "./.." "./../.."

try

    $ libtoolize --force
    $ aclocal
    $ autoreconf -f -i -Wall,no-obsolete
    $ ./configure
    $ make

For recompile the source is recommended to remove all the object files:

    $ make clean

And for the documentation:

    $ make doc

Testing
-------

If the building process have finish without error the file `PCx` must be
create. Lets try it with the smallest mps file, `afiro.mps`.

    $ ./PCX netlib-afiro

Advanced
--------

Some options MUST be set when configuring the build environment using:

    --enable-debug          Enable debug features
    --enable-warnings       Enable warnings
    --enable-cgm            Enable conjugate gradient method

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
