Presolver
=========

Linear programming models frequently contain redundant information, as well as
other information and structure that allow some components of the solution to be
determined without recourse to a sophisticated algorithm. The purpose of
presolve or prepocessing routines is to detect and handle these features of the
input, producing a (smaller) problem to be solved by the actual linear
programming algorithm. Presolvers signigicantly enhance the efficiency and
robustness of both simplex and interior-point codes.

The presolver in PCx works with the formulation stored in the ``LPtype`` data
structure. It makes use of techniques described by Andersen and Andersen,
checking the data for the following features:

Infeasibility
  Check that :math:`u_i \geq 0` for each upper bound :math:`u_i`, :math:`i \in
  \mathcal{U}`, and taht a zero of :math:`A` has a corresponding zero in the
  right-hand side vector :math:`b`.

Empty Rows
  If the matrix :math:`A` has a zero row and a corresponding zero in the
  :math:`b`, it can be removed from the problem.

Duplicate Rows
  When a row :math:`A` (and the corresponding element of the right-hand side
  :math:`b`) is simply a multiple of another row, we can delete it without
  affecting the primal solution.

Empty Columns
  The correspoding element :math:`x_i` can be fixed at either its lower or upper
  bound, depending of the sign of the cost vector coefficient :math:`c_i`. If
  the required bound does not exist, the problem is declared to be primal
  unbounded.

Fixed variables.
  If the variables has lower and upper bounds both zero, it can obviously be
  fixed at zero and removed from the problem.

Singleton Rows
  If the :math:`i`-th row of :math:`A` contains the single nonzero element
  :math:`A_{ij}`, we clearly have :math:`x_j = b_i / A_{ij}`, so this variable
  can be removed from the problem. The :math:`i`-th row of :math:`A` (and hence
  the dual variable :math:`\pi_i`) can also be removed.

Singlethon Columns
  When :math:`A_{ij}` is the only nonzero in the column :math:`j` of :math:`A`,
  and :math:`x_j` is a free variable, we can express :math:`x_j` in terms of the
  ohter variables represented in row :math:`i` of :math:`A` and eliminate it
  from the problem. Even if not free, :math:`x_j` can be eliminated if its
  bounds are weaker than those implied by the ranges of the other elements
  represented in the row :math:`A_i`.

Forced Rows
  Some times, the linear constraint represented by row :math:`i` of :math:`A`
  forces all its variables to either their upper or lower bounds. An example
  would be the constraint :math:`10 x_3 - 4 x_{10} + x_{12} = -4` subject to the
  bounds :math:`x_3 \in [ 0, +\infty )`, :math:`x_{10} \in [0, 1]`,
  :math:`x_{12} \in [0, +\infty )`. In this case, we must have :math:`x_3 = 0`,
  :math:`x_{10} = 1`, and :math:`x_{12} = 0`, so these three variables (and the
  corresponding row of :math:`A`) can be eliminated.

The presolver makes multiples passes throught the data, checking for each of the
above features in turn. Problem reductions on one pass frequently uncover futher
reductions that are detected on subsequent passes. The presolver terminates when
a complete pass is performed without detectiong futther opportunities for
reduction. Each reduction operation is pushed onto a stack, which is
subsequently popped after the solution of the reduced linear program is found.
The effect of popping the stack is to express the solution in terms of the
original, unreduced formulation.

Despite the complexity of the code, the presolver requires little CPU time in
comparison with a single iteration of the interior-point solver.

Code for the presolver can be found in the file ``presolver.c``. The data
structures are defined in ``pre.h``. This code can be used on a stand-alone
basis independently of the PCx solver, to presolver any linear program supplied
in the ``LPtype`` format.
