Data Structures
===============

PCx has it own data structures, many of then can be find in ``main.h``, and here
you will find some explanation about it.

.. important::
   The matrix used by PCx are stored in the usual compressed-sparse-row (CSR)
   format.

.. warning::
   PCx uses Fortran-style indexing, in which row and column indices start at 1
   rather than zero.

``sparseMatrix``
----------------

The ``sparseMatrix`` data structure contains a single matrix in the CSR format.

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

it will be stored as

============= === === === === === === === === === === === === ===
Pointers
============= === === === === === === === === === === === === ===
``pBeginRow``   1   4   6   9  12
``pEndRow``     4   6   9  12  14
``Row``         1   2   4   1   2   3   4   5   1   3   4   2   5
``Value``       1  -1  -3  -2   5   4   6   4  -4   2   7   8  -5
============= === === === === === === === === === === === === ===

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

``Factor Type``
---------------

.. literalinclude:: ../SRC/main.h
   :language: c
   :lines: 94-104

