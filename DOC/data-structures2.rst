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

Before the execution of ``ordmmd``:

.. literalinclude:: gdb/maros0_ordmmd_before.txt
   :end-before: +set loggin off

After the execution of ``ordmmd``:

.. literalinclude:: gdb/maros0_ordmmd_after.txt
   :end-before: +set loggin off

The sample bellow are based in ``afiro.mps`` (download it :download:`here
<_static/afiro.mps>`).

Before the execution of ``ordmmd``:

.. literalinclude:: gdb/afiro_ordmmd_before.txt
   :end-before: +set loggin off

After the execution of ``ordmmd``:

.. literalinclude:: gdb/afiro_ordmmd_after.txt
   :end-before: +set loggin off

``sfinit``
^^^^^^^^^^

The sample bellow are based in ``maros0.mps`` (download it :download:`here
<_static/maros0.mps>`).

Before the execution of ``sfinit``:

.. literalinclude:: gdb/maros0_sfinit_before.txt
   :end-before: +set loggin off

After the execution of ``sfinit``:

.. literalinclude:: gdb/maros0_sfinit_after.txt
   :end-before: +set loggin off

The sample bellow are based in ``afiro.mps`` (download it :download:`here
<_static/afiro.mps>`).

Before the execution of ``sfinit``:

.. literalinclude:: gdb/afiro_sfinit_before.txt
   :end-before: +set loggin off

After the execution of ``sfinit``:

.. literalinclude:: gdb/afiro_sfinit_after.txt
   :end-before: +set loggin off

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
