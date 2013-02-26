Invoking
========

To solve a linear program contained in the MPS file ``probname.mps``, one should
go to the working directory (that is, the directory in which the executable
``PCx`` resides) and type ::

    $ ./PCx probname

The file ``probname.mps`` can reside either in the working directory or in an
"input directory" defined in the specifications file. PCx first search the input
directory (if specified) for the given file. It search for the file name both
with and without the ``.mps`` extension, If it does not find the file in the
input directory, it search the working directory.

PCx optionally produces two output files named ``probname.out`` and
``probname.log``, according to the options supllied by the user in the
specifications file. These files are written in the working directory. They
contain, respectively, the primal-dual point returned by the algorithm (provied
the termination status is not ``infeasible``), and a summary of the iteration
history, timmings, preprocessor results, and sparsity statistics for the
Cholesky factorization. Output is also written to standard output during
execution of PCx. Essentially, the on-screen output consists of the information
written to the file ``probname.log``, together with error messages and warnings.

When PCx is executed as a standalone system and a runtime error is detected, the
code returns a nonzero integer to the operating system. The return status
indicates the type of error, as follows:

``1``
  invocation error for ``PCx``.

``2``
  memmory allocation error (usually, insufficient storage available).

``3``
  error in the MPS input file.

``4``
  error in the specifications file.

``5``
  errror detected during presolve.

``6``
  error encountered during matrix factorization, conjugate gradient iteration,
  sparse matrix multiplication, or dense column linear algebra.

The subroutine ``PCx()`` can also be invoked directly from user-written code. In
this case, the user should fill out data structures that define the linear
program and the algorihmic parameters.
