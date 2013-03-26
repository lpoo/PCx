Specifications File
===================

PCx allow many algorithmic parameters and options to be set by the user. These
quantities are stored internally via a data structure of type ``Parameters``.

If the user provides input to PCx via an MPS file (rather than involking
``PCx()`` directly via a subroutine call), the ``Parameters`` data structure is
allocated automatically by the program and default values are assigned to all
parameters. You can override the default values by defining a specifications
file, which contaions a number of keywords and numerical values.

PCx search for the specificaitions file in a number of locations. If the names
of the MPS input file is ``probname.mps``, PCx look for the following files, in
order:

#. ``probname.spc``
#. ``probname.specs``
#. ``spc``
#. ``specs``
#. ``PCx.specs``

If more than one of these file exist, PCx uses the first file in the list above
and ignores the others.

The following is a list of keywords, in alphabetic order, that can be used in
the specifications file, together with theri default settings. The file should
contain one such keyword per line, together with its correspondigin numerical
value or option, if appropriate. The file is processed sequentially from top to
bottom, so the effect of nay line in the file can be undone by a later line. For
keywords with a ``yes/no`` argument, omission of the argument will be taken to
mean ``yes``.  (The default setting is not necessarily ``yes``.) In the
descriptions below, we assume that PCx is involed with the comand ::

    $ ./PCx probname

``boundname {name}``
  Request the bound to the specific column ``name`` in ``probname.mps``.

  **Default**: the first BOUND in the MPS file is used.

``cachesize {value}``
  Input the size of the cache on the machine, in Kilobytes. Any value in the
  range 0-2048 is acceptable. Specify 0 for Cray machines. This parameters is
  used by the Ng-Peyton sparse Cholesky code.

  **Default**: 16.

``centerexp {value}``
  Specify the exponent to be used for calculation of the centering parameter.
  Any real value in the range 1.0-4.0 is allowable.

  **Default**: 3.0.

``dualfeastol {value}``
  Specify a dual feasibility tolerance.

  **Default**: 10^{-8}.

``history {yes}/{no}``
  Request that a history file be written (``yes``) or not written (``no``). If
  the ``yes``, the file ``probname.log`` is written to the working directory.

  **Default**: ``no``.

``HOCorrections {yes}/{no}``
  Request that Gondzio's higher-order corrections be used to enhance the search
  direction.

  **Default**: ``yes``.

``inputdirectory {name}``
  If PCx is to search for the MPS input file in another directory, in addition
  to the current working directory, name this other directory here. Remember to
  include a trailing ``/``. PCx always looks first in the current working
  directory. If it cannot find the file there, it looks in the specified input
  directory. The output and history files always are written to the working
  directory.

  **Defalt**: working directory.

``iterationlimit {value}``
  An upper limit on the number of iterations. Any positive integer is allowable.

  **Default**: 100.

``max``
  Maximize the objective.

  **Default**: ``no``.

``MaxCorrections {value}``
  If ``HOCorrections = yes``, the parameter ``MaxCorrections`` is an upper limit
  on the number of Gondzio's higher-order corrections allowed at each iteration.
  If ``value = 0``, the maximum is determined automatically by PCx according to
  the relative cost of factorization and solve operations. If ``HOCorrections =
  no``, ``MaxCorrections`` is ignored.

  **Default**: 0.

``mim``
  Minimize the objective.

  **Default**: ``yes``.

``objectivename {name}``
  Request the objective cost vector to be the specific row ``name`` in
  ``probname.mps``.

  **Default**: the first row of type ``N`` in ``probname.mps`` is take to be the
  objective.

``optol {value}``
  Specify an optimality tolerance.

  **Default**: 10^{-8}.

``OrderAlg {value}``
  Specify the ordering algorithm to be use. The list of algorithm avaliable are:

  0. The natural order are used.
  1. The multiple minimum degree ordering.
  2. The Reverse Cuthill-McKee.

  **Default**: 1.

``preprocess {yes}/{no}``
  Synonymous with ``presolver``.

  **Default**: ``yes``.

``presolve {yes}/{no}``
  Request that presolving be performed (``yes``) or not preformed (``no``).

  **Deafult**: ``yes``.

``prifeastol {yes}/{no}``
  Specify a primal feasibility tolerance.

  **Default**: 10^{-8}.

``rangename {name}``
  Request the range to be the specific column ``name`` in ``probname.mps``.

  **Default**: the first range encountered in the MPS file is used.

``refinment {yes}/{no}``
  Perform preconditioned conjugate gradient refinement of the computed solution
  to the linear system if it has a relative residual larger than the parameter
  ``prifeastol`` (``yes``) or don't perform any iterative refinement (``no``).

  **Default**: ``no``.

``rhsname {name}``
  Request the right-hand side to be the specific column ``name`` in
  ``probname.mps``.

  **Default**: the first RHS encountered in the MPS file is used.

``scaling {yes}/{no}``
  If ``yes``, row and column scaling is performed on the constraint matrix.

  **Default**: ``yes``.

``solution {yes}/{no}``
  Request that a solution file be written (``yes``) or not written (``no``). If
  the solution file is written, it is named ``probname.out`` and is placed in
  the working directory.

  **Default**: ``yes``.

``stepfactor {value}``
  Specify a value in the range (0, 1) that is usde in Mehrotra's adaptive
  steplength heuristic. This value is a lower bound for $gamma^P$ and $gamma^D$.

  **Default**: 0.9.

``unrollinglevel {value}``
  Specify the level of loop unrolling. Allowable values are 1, 2, 4, and 8.
  (This parameter is used only in the Ng-Peyton sparse Cholesky code).

  **Default**: 4.

If you call ``PCx()`` directly from your own code, you must fill out the
``Parameters`` data structure explicitly. This task is easier if you use the
routine ``*NewParameters()`` to allocate the storage, since this routine assigns
default values to all the parameters. You can then make any desired alterations
before passing the data structure to the ``PCx()`` routine.
