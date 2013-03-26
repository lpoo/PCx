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

To test PCx it on one of the input files int the directory ``mps``, modify the
sample specifications file ``PCx.specs`` if desired, the type ::

    $ ./PCX afiro

or ::

    $ ./PCx 25fc47
