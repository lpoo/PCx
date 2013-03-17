Data Structures (2)
===================

In this section you will find some explanation about the data structures used in
the Ng-Peyton routines for solving sparse symmetric positive definite linear
systems on a single processor using sparse Cholesky factorization.

Fortran Implementation
----------------------

The Ng-Peyton routines was written in Fotran and can be found in the directory
``Ng-Peyton``. Below you will find the explanation about some variables names
used in Ng-Peyton.

``NEQNS``
    There are ``NEQNS`` columns (rows) in the matrix :math:`A` (also :math:`L`).

``INVP(1:NEQNS)``
    Contains the inverse permutation vector, which maps old locations to new
    locations.  More specifically, ``INVP(J)`` is the new location of column
    ``J`` (row ``J``) of a after the reordering has been applied to the columns
    and rows of the matrix.

``NSUPER``
    There are ``NSUPER`` supernodes.  The supernodes are numbered ``1, 2, ... ,
    NSUPER``.  It is important to note that ``NSUPER <= NEQNS``.

``SNODE(1:NEQNS)``
    ``SNODE(J)`` is the supernode to which column ``J`` OF :math:`L` belongs.

``XSUPER(1:NSUPER+1)``
    For ``1 <= JSUP <= NSUPER``, ``XSUPER(JSUP)`` is the first column of
    supernode ``JSUP``.  Supernode ``JSUP`` contains columns ``XSUPER(JSUP),
    XSUPER(JSUP)+1, ... , XSUPER(JSUP+1)-1``.

``NSUB``
    Number of row subscripts stored in ``LINDX(*)``.

``LINDX(1:NSUB)``
    ``LINDX(*)`` contains, in column major order, the row subscripts 
    of the nonzero entries in :math:`L` in a compressed storage format.
    This format is described immediately below.

``XLINDX(1:NSUPER)``
    Consider supernode ``JSUP``, and let ``FSTCOL = XSUPER(JSUP)`` - the 
    first column in supernode ``JSUP``.  Let ``JCOL`` be any column in
    supernode ``JSUP`` (possibly ``JCOL = FSTCOL``).  Let 
    ``FSTSUB = XLINDX(JSUP)+JCOL-FSTCOL``, and 
    ``LSTSUB = XLINDX(JSUP+1) - 1``.  The row subscripts of the 
    nonzero entries in column ``JCOL`` of :math:`L` lie in ascending order 
    in ``LINDX``; they are ``LINDX(FSTSUB), LINDX(FSTSUB+1), ... , 
    LINDX(LSTSUB)``.  Note that a row subscript for the main 
    diagonal entry is included, i.e., ``LINDX(FSTSUB) = JCOL``.

``NNZL``
    Number of nonzero entries in :math:`L`.

``LNZ(1:NNZL)``
    ``LNZ(*)`` contains, in column major order, the nonzero entries
    in :math:`L`.

``XLNZ(1,NEQNS+1)``
    For ``1 <= J <= NEQNS``, ``XLNZ(J)`` points to the first nonzero 
    entry in column ``J`` of :math:`L`.  The nonzero entries in column ``J`` OF :math:`L` 
    lie in ascending order by row subscript in ``LNZ(*)``; they are 
    ``LNZ(XLNZ(J)), LNZ(XLNZ(J)+1), ... , LNZ(XLNZ(J+1)-1)``.

C Interface
-----------

The C interface to the Ng-Peyton use the ``NGPeytonType`` data structure:

.. literalinclude:: ../SRC/Ng-Peyton.h
   :language: c
   :lines: 2-11

The following Ng-Peyton subroutines are present in the C interface: 

``ordmmd``
^^^^^^^^^^
The sample bellow are based in ``maros0.mps`` (download it :download:`here
<_static/maros0.mps>`).

Before the execution of ``ordmmd``::

    dimension                  = 3
    NgPeyton.pSuperNodeCols    = 1
    NgPeyton.SuperNodeRows     = 2
    Factor.Perm                = {0, 0, 0}
    Factor.InvPerm             = {0, 0, 0}
    NgPeyton.NumCompressedCols = 0

After the execution of ``ordmmd``::

    dimension                  = 3
    NgPeyton.pSuperNodeCols    = 1
    NgPeyton.SuperNodeRows     = 2
    Factor.Perm                = {0, 0, 0}
    Factor.InvPerm             = {0, 0, 0}
    NgPeyton.NumCompressedCols = 0

The sample bellow are based in ``afiro.mps`` (download it :download:`here
<_static/afiro.mps>`).

Before the execution of ``ordmmd``::

    dimension                  = 27
    NgPeyton.pSuperNodeCols    = 1
    NgPeyton.SuperNodeRows     = 2
    Factor.Perm                = {0 <repeats 27 times>}
    Factor.InvPerm             = {0 <repeats 27 times>}
    NgPeyton.NumCompressedCols = 30

After the execution of ``ordmmd``::

    dimension                  = 27
    NgPeyton.pSuperNodeCols    = 1
    NgPeyton.SuperNodeRows     = 2
    Factor.Perm                = {0 <repeats 27 times>}
    Factor.InvPerm             = {0 <repeats 27 times>}
    NgPeyton.NumCompressedCols = 30

``sfinit``
^^^^^^^^^^

The sample bellow are based in ``maros0.mps`` (download it :download:`here
<_static/maros0.mps>`).

Before the execution of ``sfinit``::

    dimension                     = 3
    nonzeros                      = 6
    TempBeginRow                  = {1, 3, 5}
    TempRow                       = {2, 3, 1, 3, 1}
    Factor.Perm                   = {3, 1, 2}
    Factor.InvPerm                = {2, 3, 1}
    ColumnnCount                  = {0, 0, 0}
    Factor.NonzerosL              = 0
    NgPeyton.NumCompressedCols    = 2
    NgPeyton.NumSuperNodes        = 0
    NgPeyton.mapColumnToSupernode = {0, 0, 0}
    NgPeyton.SuperPartitioning    = {0, 0, 0}

After the execution of ``sfinit``::

    dimension                     = 3
    nonzeros                      = 6
    TempBeginRow                  = {1, 3, 5}
    TempRow                       = {2, 3, 1, 3, 1}
    Factor.Perm                   = {3, 1, 2}
    Factor.InvPerm                = {2, 3, 1}
    ColumnnCount                  = {3, 2, 1}
    Factor.NonzerosL              = 6
    NgPeyton.NumCompressedCols    = 3
    NgPeyton.NumSuperNodes        = 1
    NgPeyton.mapColumnToSupernode = {1, 1, 1}
    NgPeyton.SuperPartitioning    = {1, 4, 2}

The sample bellow are based in ``afiro.mps`` (download it :download:`here
<_static/afiro.mps>`).

Before the execution of ``sfinit``::

    dimension                     = 27
    nonzeros                      = 126
    TempBeginRow                  = {1, 6, 10, 13, 15, 23, 30, 34, 38, 42, 46, 52, 56, 59, 61, 68, 76, 80, 84, 88, 92, 101, 105, 112, 116, 123, 125}
    TempRow                       = {2, 3, 4, 22, 24, 1, 3, 24, 26, 1, 2, 24, 1, 5, 4, 6, 7, 8, 9, 10, 23, 25, 5, 7, 8, 9, 10, 25, 27, 5, 6, 21, 25, 5, 6, 21, 25, 5, 6, 21, 25, 5, 6, 21, 25, 12, 13, 14, 21, 22, 24, 11, 13, 22, 26, 11, 12, 22, 11, 16, 16, 17, 18, 19, 20, 23, 27, 14, 15, 17, 18, 19, 20, 23, 25, 15, 16, 21, 23, 15, 16, 21, 23, 15, 16, 21, 23, 15, 16, 21, 23, 7, 8, 9, 10, 11, 17, 18, 19, 20, 1, 11, 12, 13, 5, 15, 16, 17, 18, 19, 20, 1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 16, 2, 12, 6}
    Factor.Perm                   = {27, 26, 14, 4, 13, 3, 12, 24, 22, 2, 1, 11, 20, 19, 18, 17, 10, 9, 8, 7, 15, 25, 23, 6, 21, 5, 16}
    Factor.InvPerm                = {11, 10, 6, 4, 26, 24, 20, 19, 18, 17, 12, 7, 5, 3, 21, 27, 16, 15, 14, 13, 25, 9, 23, 8, 22, 2, 1}
    ColumnnCount                  = {0 <repeats 27 times>}
    Factor.NonzerosL              = 0
    NgPeyton.NumCompressedCols    = 79
    NgPeyton.NumSuperNodes        = -147351832
    NgPeyton.mapColumnToSupernode = {0 <repeats 27 times>}
    NgPeyton.SuperPartitioning    = {0 <repeats 27 times>}

After the execution of ``sfinit``::

    dimension                     = 27
    nonzeros                      = 126
    TempBeginRow                  = {1, 6, 10, 13, 15, 23, 30, 34, 38, 42, 46, 52, 56, 59, 61, 68, 76, 80, 84, 88, 92, 101, 105, 112, 116, 123, 125}
    TempRow                       = {2, 3, 4, 22, 24, 1, 3, 24, 26, 1, 2, 24, 1, 5, 4, 6, 7, 8, 9, 10, 23, 25, 5, 7, 8, 9, 10, 25, 27, 5, 6, 21, 25, 5, 6, 21, 25, 5, 6, 21, 25, 5, 6, 21, 25, 12, 13, 14, 21, 22, 24, 11, 13, 22, 26, 11, 12, 22, 11, 16, 16, 17, 18, 19, 20, 23, 27, 14, 15, 17, 18, 19, 20, 23, 25, 15, 16, 21, 23, 15, 16, 21, 23, 15, 16, 21, 23, 15, 16, 21, 23, 7, 8, 9, 10, 11, 17, 18, 19, 20, 1, 11, 12, 13, 5, 15, 16, 17, 18, 19, 20, 1, 2, 3, 11, 5, 6, 7, 8, 9, 10, 16, 2, 12, 6}
    Factor.Perm                   = {14, 4, 3, 24, 13, 26, 12, 22, 2, 1, 11, 10, 9, 8, 7, 25, 20, 19, 18, 17, 27, 15, 23, 6, 21, 5, 16}
    Factor.InvPerm                = {10, 9, 3, 2, 26, 24, 15, 14, 13, 12, 11, 7, 5, 1, 22, 27, 20, 19, 18, 17, 25, 8, 23, 4, 16, 6, 21}
    ColumnnCount                  = {3, 3, 4, 4, 4, 3, 4, 4, 3, 3, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 4, 3, 2, 1}
    Factor.NonzerosL              = 107
    NgPeyton.NumCompressedCols    = 94
    NgPeyton.NumSuperNodes        = 22
    NgPeyton.mapColumnToSupernode = {1, 2, 3, 4, 5, 6, 7, 8, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 22, 22, 22, 22}
    NgPeyton.SuperPartitioning    = {1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 28, 5, 4, 3, 2}

``bfinit``
^^^^^^^^^^

``blklvl``
^^^^^^^^^^

``blkslv``
^^^^^^^^^^

``blkslf``
^^^^^^^^^^

``blkslb``
^^^^^^^^^^

``symfct``
^^^^^^^^^^

``inpnv``
^^^^^^^^^

``ordnat``
^^^^^^^^^^
