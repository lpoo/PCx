Data Structures
===============

PCx has it own data structures, many of then can be find in ``main.h``, and here
you will find some explanation about it.

.. important::
   The matrix used by PCx are stored in the usual compressed-sparse-row (CSR)
   format. For more information about sparse matrix storage see the references.

.. warning::
   PCx uses Fortran-style indexing, in which row and column indices start at 1
   rather than zero.

``MMTtype``
-----------

The ``MMTtype`` data structure contains a sparse structure with pointer
information for both the matrix and its transpose.

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 73-80


``sparseMatrix``
----------------

The ``sparseMatrix`` data structure contains a single matrix in the CSR format
(for more information see [Intel]_).

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 84-90

Given the matrix B,

.. math::

   B &= \begin{bmatrix}
        1 & -1 & 0 & -3 & 0 \\
        -2 & 5 & 0 & 0 & 0 \\
        0 & 0 & 4 & 6 & 4 \\
        -4 & 0 & 2 & 7 & 0 \\
        0 & 8 & 0 & 0 & -5
        \end{bmatrix},

it will be stored as ::

    NumRows = 5,
    NumCols = 5,
    Nonzeros = 13,
    pBeginRow = {1, 4, 6, 9, 12},
    pEndRow = {4, 6, 9, 12, 14},
    Row = {1, 2, 4, 1, 2, 3, 4, 5, 1, 3, 4, 2, 5},
    Value = {1, -1, -3, -2, 5, 4, 6, 4, -4, 2, 7, 8, -5}

``Factor Type``
---------------

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 94-104

Sample
^^^^^^

The sample bellow are based in ``maros0.mps`` (download it :download:`here
<_static/maros0.mps>`). ::

    ptr = 0x62da00,
    FactorizationCode = 0x62d810 "Ng-Peyton sparse Cholesky library",
    N = 3,
    NonzerosL = 6,
    SmallDiagonals = 0,
    Perm = {3, 1, 2},
    InvPerm = {2, 3, 1},
    AAT = {
        NumRows = 3,
        NumCols = 3,
        Nonzeros = 9,
        pBeginRow = {1, 4, 7},
        pBeginRowT = 0x0,
        pEndRow = {3, 6, 9},
        pEndRowT = 0x0,
        Row = {1, 2, 3, 1, 2, 3, 1, 2, 3},
        RowT = 0x0,
        Value = {0, 0, 0, 0, 0, 0, 0, 0, 0},
        ValueT = 0x0
    },
    Ndense = 0,
    maskDense = 0x62d640,
    W = 0x62d8e0,
    Ldense = 0x62d9c0

``MPStype``
-----------

The ``MPStype`` data structure contains a complete specification of a single
linear programming problem in the general formulation, i.e., may include upper
and lower bounds, linear equality constraints, linear inequality constraints and
free variables. This data strucutre also stores the names assigned to the rows,
columns, and objectives of the model specified in the MPS file.

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 108-140

``MPSchanges``
--------------

``LPtype``
----------

The ``LPtype`` data structure contains a single linear program in the follow
form:

.. math::

   \min_{x in \mathbb{R}^n} c^T x \text{ subject to } A x = b,
                                                      0 \leq x_i, i \in \mathcal{N},
                                                      0 \leq x_i \leq u_i, i \in \mathcal{U},
                                                      x_i \text{ free}, i \in \mathcal{F}.

Also, the ``LPtype`` data structure are used to store linear program in the form

.. math::

   \min_{x in \mathbb{R}^n} c^T x \text{ subject to } A x = b,
                                                      0 \leq x_i, i \in \mathcal{N},
                                                      0 \leq x_i \leq u_i, i \in \mathcal{U},

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 161-194

Sample
^^^^^^

The sample bellow are based in ``maros0.mps`` (download it :download:`here
<_static/maros0.mps>`). :: 

    Rows = 3,
    Cols = 7,
    Ents = 13,
    Atranspose = {
        NumRows = 7,
        NumCols = 3,
        Nonzeros = 13,
        pBeginRow = {1, 5, 9, 0, 0, 0, 33},
        pEndRow = {4, 8, 13, 0, 0, 0, 33},
        Row = {1, 3, 4, 5, 2, 3, 4, 6, 1, 2, 3, 4, 7},
        Value = {1, 1, 1.5, 1, 1.5, 0.5, 0.5, 1, 2.5, 2, 3, 2, -1}
    },
    A = {
        NumRows = 3,
        NumCols = 7,
        Nonzeros = 13,
        pBeginRow = {1, 3, 5, 8, 11, 12, 13},
        pEndRow = {2, 4, 7, 10, 11, 12, 13},
        Row = {1, 3, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3},
        Value = {1, 2.5, 1.5, 2, 1, 0.5, 3, 1.5, 0.5, 2, 1, 1, -1}
    },
    b = {50, 35, 125},
    c = {4.5, 2.5, 4, 4, 0, 0, 0},
    cshift = -40,
    VarType = {0, 0, 2, 2, 0, 0, 2},
    UpBound = {0, 0, 30, 25, 0, 0, 10},
    NumberBounds = 2,
    BoundIndex = 0x62a080,
    NumScale = 0,
    ColScale = 0x0,
    RowScale = 0x0,
    NumberFree = 0,
    FreeIndex = 0x62a0b0,
    NumberSplit = 0,
    FreePlus = 0x0,
    FreeMinus = 0x0

``Parameters``
--------------

``IterationRecord``
-------------------

``FactorizationRecord``
-----------------------

``solution``
------------

``Iterate``
-----------

``GenralInfo``
--------------

.. rubric:: Refereces

.. [Intel] Intel Corporation. Sparse Matrix Storage Formats. http://software.intel.com/sites/products/documentation/hpc/mkl/mklman/GUID-9FCEB1C4-670D-4738-81D2-F378013412B0.htm.

.. [Netlib] Jack Dongarra. Survey of Sparse Matrix Storage Formats. http://netlib.org/utk/papers/templates/node90.html.
