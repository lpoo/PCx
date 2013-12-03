Formulation
===========

PCx accepts any valid linear program that can be specified in the MPS format.

.. note::

   This description of the MPS file format have been copied from
   http://web.mit.edu/lpsolve/doc/mps-format.htm

Fixed MPS format
----------------

The main things to know about fixed MPS format are that it is column oriented
(as opposed to entering the model as equations), and everything
(variables, rows, etc.) gets a name.

MPS is an old format, so it is set up as though you were using punch
cards. Fields start in column 2, 5, 15, 25, 40 and 50.
Sections of an MPS file are marked by so-called header cards,
which are distinguished by their starting in column 1. Although it is
typical to use upper-case throughout the file (like I said, MPS has
long historical roots), many MPS-readers will accept mixed-case for
anything except the header cards, and some allow mixed-case anywhere.
The names that you choose for the individual entities (constraints or
variables) are not important to the solver; you should pick names that
are meaningful to you, or will be easy for a post-processing code to
read.

Here is a little sample model written in MPS format (explained in more
detail below)::

    NAME          TESTPROB
    ROWS
     N  COST
     L  LIM1
     G  LIM2
     E  MYEQN
    COLUMNS
        XONE      COST                 1   LIM1                 1
        XONE      LIM2                 1
        YTWO      COST                 4   LIM1                 1
        YTWO      MYEQN               -1
        ZTHREE    COST                 9   LIM2                 1
        ZTHREE    MYEQN                1
    RHS
        RHS1      LIM1                 5   LIM2                10
        RHS1      MYEQN                7
    BOUNDS
     UP BND1      XONE                 4
     LO BND1      YTWO                -1
     UP BND1      YTWO                 1
    ENDATA

For comparison, here is the same model written out in lp-format::

    min: +XONE +4 YTWO +9 ZTHREE;
    LIM1: +XONE +YTWO <= 5;
    LIM2: +XONE +ZTHREE >= 10;
    MYEQN: -YTWO +ZTHREE = 7;
    XONE <= 4;
    YTWO >= -1;
    YTWO <= 1;

Strangely, there is nothing in MPS format that specifies the direction
of optimisation.  And there really is no standard "default" direction;
some LP codes will maximize if you don't specify otherwise, others will
minimize, and still others put safety first and have no default and
require you to specify it somewhere in a control program or by a
calling parameter. lp_solve uses minimization as default. This can be
changed do maximization by calling set_maxim or set_sense after reading
the file or by using option -max on the lp_solve command line program.
If you have a model formulated for minimization and the code you are
using insists on maximization (or vice versa), it may be easy to convert:
just multiply all the coefficients in your objective function by (-1).
The optimal value of the objective function will then be the negative of
the true value, but the values of the variables themselves will be correct.

Any line with an asterisk (*) in Column 1 is treated as a comment.

The eight character names used to specify variables, constraints and
other entities are fixed format. Names are not automatically justified,
so blanks are treated just like other characters. For example "ROW1    "
is not the same as " ROW1   ". (Note that some optimisers do not permit
blanks in names.) No case conversion is performed, so "row1    " is
different from "ROW1    ".

Floating point numbers may be specified in free format within the 12
character field (including embedded blanks). The following list describes
the possible ways of writing a number.

Mantissa:

``+`` or ``-``
    optional sign character (no sign indicates a positive number)
digits
    optional integer part of the mantissa
``.``
    optional decimal point (if not present, a decimal point will be assumed
    after the mantissa digit)
digits
    optional fraction part of the mantissa -the mantissa must contain at least
    one digit

Exponent (optional):

``D`` or ``E``
    exponent leader
``+`` or ``-``
    optional exponent sign
digits
    exponent digits

Numbers with an absolute value greater than 1010 or less than 10-10 are rejected.

The NAME card can have anything you want, starting in column 15.  The
ROWS section defines the names of all the constraints; entries in
column 2 or 3 are E for equality rows, L for less-than ( <= ) rows, G
for greater-than ( >= ) rows, and N for non-constraining rows (the
first of which would be interpreted as the objective function).  The
order of the rows named in this section is unimportant.

The largest part of the file is in the COLUMNS section, which is the
place where the entries of the A-matrix are put. All entries for a given
column must be placed consecutively, although within a column the
order of the entries (rows) is irrelevant. Rows not mentioned for a
column are implied to have a coefficient of zero.

The RHS section allows one or more right-hand-side vectors to be
defined; most people don't bother having more than one.  In the above
example, the name of the RHS vector is RHS1, and has non-zero values
in all 3 of the constraint rows of the problem.  Rows not mentioned in
an RHS vector would be assumed to have a right-hand-side of zero.
Note that the objective may also have a constant. This can also be
specified in this section by using the object name as constraint name
and then specifying the constant. Note that there are 2 interpretations
of this constant. Some solvers see this as the constant that would be
really in the RHS and when brought into the objective (LHS), it is negated.
Other solvers, as lp_solve does, use the specified value in the MPS file
as the value for the objective and don't negate it.

The optional BOUNDS section lets you put lower and upper bounds on
individual variables (no * wild cards, unfortunately), instead of
having to define extra rows in the matrix.  All the bounds that have
a given name in column 5 are taken together as a set.  Variables not
mentioned in a given BOUNDS set are taken to be non-negative (lower
bound zero, no upper bound).  A bound of type UP means an upper bound
is applied to the variable.  A bound of type LO means a lower bound is
applied.  A bound type of FX ("fixed") means that the variable has
upper and lower bounds equal to a single value.  A bound type of FR
("free") means the variable has neither lower nor upper bounds.

There is another optional section called RANGES that I won't go into
here. The final card must be ENDATA, and yes, it is spelled funny.

MPS input format was originally introduced by IBM to express linear
and integer programs in a standard way.  The format is a fixed column
format, so care must be taken that all information is placed in the
correct columns as described below.

The following is not intended as a complete description of MPS format,
but only as a brief introduction.  For more information, the reader is
directed to:

* "Advanced Linear Programming," by Bruce A. Murtagh
* "Computer Solutions of Linear Programs," by J.L. Nazareth

It may be useful to look at an example MPS file while reading this
MPS information.

The following template is a guide for the use of MPS format::

    Field:    1           2          3         4         5         6
    Columns:  2-3        5-12      15-22     25-36     40-47     50-61

              NAME   problem name

              ROWS

               type     name

              COLUMNS
                       column       row       value     row      value
                        name        name                name
              RHS
                        rhs         row       value     row      value
                        name        name                name
              RANGES
                        range       row       value     row      value
                        name        name                name
              BOUNDS

               type     bound       column    value
                        name        name

              SOS
               type     CaseName    SOSName   SOSpriority
                        CaseName    VarName1  VarWeight1
                        CaseName    VarName2  VarWeight2

                        CaseName    VarNameN  VarWeightN

              ENDATA

NOTES:

A. In the ROWS section, each row of the constraint matrix must have a
   row type and a row name specified.  The code for indicating row type
   is as follows:

+------+----------------------+
| type |    meaning           |
+------+----------------------+
| E    |equality              |
+------+----------------------+
| L    |less than or equal    |
+------+----------------------+
| G    |greater than or equal |
+------+----------------------+
| N    |objective             |
+------+----------------------+
| N    |no restriction        |
+------+----------------------+

B. In the COLUMNS section, the names of the variables are defined along
   with the coefficients of the objective and all the nonzero constraint
   matrix elements.  It is not necessary to specify columns for slack or
   surplus variables as this is taken care of automatically.

C. The RHS section contains information for the right-hand side of the problem.

D. The RANGES section is for constraints of the form:  h <= constraint <= u .
   The range of the constraint is  r = u - h .  The value of r is specified
   in the RANGES section, and the value of u or h is specified in the RHS
   section.  If b is the value entered in the RHS section, and r is the
   value entered in the RANGES section, then u and h are thus defined:

+----------+-----------+---------+---------+
| row type | sign of r | h       |   u     |
+----------+-----------+---------+---------+
| G        | + or -    | b       | b + |r| |
+----------+-----------+---------+---------+
| L        | + or -    | b - |r| | b       |
+----------+-----------+---------+---------+
| E        | +         | b       | b + |r| |
+----------+-----------+---------+---------+
| E        | -         | b - |r| | b       |
+----------+-----------+---------+---------+

E. In the BOUNDS section, bounds on the variables are specified.  When
   bounds are not indicated, the default bounds ( 0 <= x < infinity )
   are assumed.  The code for indicating bound type is as follows:

+-----+--------------------------------------------+
|type |          meaning                           |
+-----+--------------------------------------------+
| LO  | lower bound        b <= x (< +inf)         |
+-----+--------------------------------------------+
| UP  | upper bound        (0 <=) x <= b           |
+-----+--------------------------------------------+
| FX  | fixed variable     x = b                   |
+-----+--------------------------------------------+
| FR  | free variable      -inf < x < +inf         |
+-----+--------------------------------------------+
| MI  | lower bound -inf   -inf < x (<= 0)         |
+-----+--------------------------------------------+
| PL  | upper bound +inf   (0 <=) x < +inf         |
+-----+--------------------------------------------+
| BV  | binary variable    x = 0 or 1              |
+-----+--------------------------------------------+
| LI  | integer variable   b <= x (< +inf)         |
+-----+--------------------------------------------+
| UI  | integer variable   (0 <=) x <= b           |
+-----+--------------------------------------------+
| SC  | semi-cont variable x = 0 or l <= x <= b    |
|     | l is the lower bound on the variable       |
|     | If none set then defaults to 1             |
+-----+--------------------------------------------+

F. Sections RANGES and BOUNDS are optional as are the fields 5 and 6.
   Everything else is required.  In regards to fields 5 and 6, consider
   the following 2 constraints:

| const1:  2x + 3y <= 6
| const2:  5x + 8y <= 20

   Two ways to enter the variable x in the COLUMNS section are::

      (Field:  2    3           4            5         6  )
    1.         x  const1       2.0         const2     5.0
    2.         x  const1       2.0
               x  const2       5.0

G. A mixed integer program requires the specification of which variables
   are required to be integer.  Markers are used to indicate the start
   and end of a group of integer variables.  The start marker has its
   name in field 2, 'MARKER' in field 3, and 'INTORG' in field 5.  The
   end marker has its name in field 2, 'MARKER' in field 3, and 'INTEND'
   in field 5.  These markers are placed in the COLUMNS section.
   When there are BOUNDS on the variables, then these are used as lower
   and upper bound of these integer variables and there is no confusion
   possible. Even a lower bound of 0 is already enough. In that case, if
   there is no upper bound, infinite is used.
   However there is an interpretation problem if there are no bounds at
   all on these variables. Some solvers then use 0 as lower bound and 1
   as upper bound. So the variables are treated as binary variables.
   That is the original IBM interpretation.
   Other solvers, like lp_solve, use the default bounds on variables in
   that case. That is 0 as lower bound and infinite as upper bound.
   When lp_solve writes an MPS file, it will write the default lower
   bound of 0 if there are no lower/upper bounds set on the variable. As
   such, there is no confusion.
   However when lp_solve reads an MPS file and there are no bounds on
   variables between INTORG/INTEND, it interpretes the variables as
   integer and not binary as some other solvers do. That could result
   in another solution than expected.

H. A specially ordered set of degree N is a collection of variables where
   at most N variables may be non-zero.  The non-zero variables must be
   contiguous (neighbours) sorted by the ascending value of their respective
   unique weights.  In lp_solve, specially ordered sets may be of any
   cardinal type 1, 2, and higher, and may be overlapping.  The number of
   variables in the set must be equal to, or exceed the cardinal SOS order.

   Below is a representation of a SOS in an MPS file, where each SOS is
   defined in its own SOS section, which should follow the BOUNDS section.  ::

    0        1         2         3         4
    1234567890123456789012345678901234567890
    SOS
     Sx CaseName  SOSName.  SOSpriority.
        CaseName  VarName1  VarWeight1..
        CaseName  VarName2  VarWeight2..
        CaseName  VarNameN  VarWeightN..

   x at the second line, position 3, defines is the order of the SOS.
   Due to limitations in the MPS format, N is restricted to the 1..9 range.
   Each SOS should be given a unique name, SOSName. lp_solve does not
   currently use case names for SOS'es and the CaseName could be any non-empty
   value.  The SOSpriority value determines the order in which multiple SOS'es
   are analysed in lp_solve.
   See also Interpolation with GAMS.
   Example::

    NAME          SOS2test
    ROWS
     N  obj
     L  c1
     L  c2
     E  c3
    COLUMNS
        x1        obj                 -1   c1                  -1
        x1        c2                   1
        x2        obj                 -2   c1                   1
        x2        c2                  -3   c3                   1
        x3        obj                 -3   c1                   1
        x3        c2                   1
        x4        obj                 -1   c1                  10
        x4        c3                -3.5
        x5        obj                  0
    RHS
        rhs       c1                  30   c2                  30
    BOUNDS
     UP BOUND     x1                  40
     LI BOUND     x4                   2
     UI BOUND     x4                   3
    SOS
     S2 SET       SOS2                10
        SET       x1               10000
        SET       x2               20000
        SET       x4               40000
        SET       x5               50000
    ENDATA

Free MPS format
---------------

The free format is very similar to the fixed MPS format, but it is less
restrictive e.g. it allows longer names. Also some implementations allow
more than 12 positions to specify the values. The fields do not have
fixed column positions as in the fixed MPS format. They may be written
anywhere except column 1, with each field separated from the next by one
or more blanks. However, they must appear in the same sequence as in
fixed format. In the rows and bounds sections, the codes can be lower
and upper case and at any starting position. Repeated column names are
sometimes skipped and spaces are put there instead. The Fortran D
exponent is allowed in the values.

There is one important limitation compared to the fixed MPS format:
names may not contain blanks.

Note that the free MPS parser cannot read all fixed MPS formats
correctly. Spaces in the names or names starting with spaces will give
problems. It is not sure that an error will be given in that case. If
the format complies to the free MPS format then it won't... So if you
know that a model is in fixed MPS format, use it and not the free
format. It is advised to first try the fixed format and only if it
doesn't work, use the free format.

Also note that there is no real standard for the free format. Each
implementer has its own implementation and interpretation of the free
format ... Some allow one space in the name, some require that names
must take at least 8 positions (and thus extended with spaces). Some
allow to have more than 6 fields on a line. lp\_solve only reads the
first 6 fields and ignores the rest.

lp\_solve tries to handle all possible free formats. The only real
limitation is that there may be no blanks in names (also no leading
blanks) and only 6 fields per line may be used. When lp\_solve writes an
mps file in free format, it will be the same as for fixed format, except
if names are longer than 8 characters. In that case all data is shifted
to the right.

OBJSENSE
^^^^^^^^

Several solvers have added a 'standard' to the free MPS format to
allow to specify the objective direction. This via the new optional
section OBJSENSE. Below this section, there may be one line that
specifies the objective direction. This in field 1 of this line via the
following possible keywords: MAXIMIZE, MAX, MINIMIZE, MIN. If the
section is not specified, then lp\_solve assumes minimization, just like
the fixed MPS format.
For example::

    OBJSENSE
     MAX

This section should be before the ROWS section.
For example::

    NAME          TESTPROB
    OBJSENSE
     MAX
    ROWS
     N  COST
     L  LIM1
     G  LIM2
     E  MYEQN
    COLUMNS
        XONE      COST                 1   LIM1                 1
        XONE      LIM2                 1
        YTWO      COST                 4   LIM1                 1
        YTWO      MYEQN               -1
        ZTHREE    COST                 9   LIM2                 1
        ZTHREE    MYEQN                1
    RHS
        RHS1      LIM1                 5   LIM2                10
        RHS1      MYEQN                7
    BOUNDS
     UP BND1      XONE                 4
     LO BND1      YTWO                -1
     UP BND1      YTWO                 1
    ENDATA

The lp\_solve free MPS reader recognises and interprets all possible
OBJSENSE direction values. When a free MPS file is created, the OBJSENSE
section will only be written when the direction is maximization. This
because minimization is by default assumed and to stay as compatible as
possible with other solvers.

OBJNAME
^^^^^^^

Several solvers have added a 'standard' to the free MPS format to
allow to specify the objective row. By default the first "N" row defined
in the ROWS section becomes a problem's objective; a different objective
may be specified in the optional OBJNAME section, which contains exactly
one data line that names the objective in field 1.
For example::

    OBJNAME
     obj2

This section should be before the ROWS section.
For example::

    NAME          TESTPROB
    OBJNAME
     PROFIT
    ROWS
     N  COST
     N  PROFIT
     L  LIM1
     G  LIM2
     E  MYEQN
    COLUMNS
        XONE      COST                -1   PROFIT               1
        XONE      LIM1                 1   LIM2                 1
        YTWO      COST                -4   PROFIT               4
        YTWO      LIM1                 1   MYEQN               -1
        ZTHREE    COST                -9   PROFIT               9
        ZTHREE    LIM2                 1   MYEQN                1
    RHS
        RHS1      LIM1                 5   LIM2                10
        RHS1      MYEQN                7
    BOUNDS
     UP BND1      XONE                 4
     LO BND1      YTWO                -1
     UP BND1      YTWO                 1
    ENDATA

The lp\_solve free MPS reader recognises and interprets this OBJNAME
section and uses the objective name specified here. Other "N" cards in
the ROWS section are then ignored. Note that if there is no OBJNAME
section that, just like in the fixed MPS format, the first "N" card from
the rows section is then taken and all other "N" cards are ignored. When
a free MPS file is created, the OBJNAME section will never be created
since lp\_solve always only has one objective function in memory.
