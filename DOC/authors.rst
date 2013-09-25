Authors
=======

The original authors of PCx are

* Joseph Czyzyk, Argonne and Poland;
* `Sanjay Mehrotra
  <http://www.iems.northwestern.edu/content/Member.asp?MemberIK=37>`_, Northwestern;
* Michael Wagner, Cornell;
* `Stephen Wright <http://www.cs.wisc.edu/~swright/>`_, Wisconsin.

PCx uses the sparse Cholesky linear algebra routines of `Esmond Ng
<http://crd-legacy.lbl.gov/~EGNg/>`_ and Barry Peyton. 

Marc Wenzel programmed the dense-column-handling and conjugate gradient
refinement features that were added for the beta-2.0 release of PCx.

Doug Moore gave valuable advice and in particular pointed out and repaired many
memory leaks in the beta-1.0 release.

Hans Mittelmann prepared the executables for numerous architectures and ran many
of the benchmark tests.

Jean-Pierre Goux set up the request form and download system. 

How to cite
-----------

To cite PCx in publications, please use:

    Czyzyk, Joseph et. al. PCx User Guide (Version 1.1). 1997.

A BibTeX entry for LaTeX users is ::

    @techreport{PCx,
        author    = {Joseph Czyzyk and Sanjay Mehrotra and Michael Wagner and
        Stephen J. Wright},
        title     = {PCx User Guide (Version 1.1)},
        institution = {Optimization Technology Center},
        year      = {1997},
    }

or ::

    @article{Czyzyk:1999hk,
        author = {Czyzyk, Joseph and Mehrotra, Sanjay and Wagner, Michael and
            Wright, Stephen J},
        title = {{PCx: an interior-point code for linear programming}},
        journal = {Optimization Methods and Software},
        year = {1999},
        volume = {11},
        number = {1},
        pages = {397--430}
    }
