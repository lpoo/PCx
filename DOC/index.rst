.. PCx documentation master file, created by
   sphinx-quickstart on Mon Feb 18 15:01:15 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PCx's unofficial documentation!
==========================================

.. note::

   A `parallel implementation of PCx
   <http://www.cs.cornell.edu/Info/People/mwagner/pPCx/>`_ was prepared by
   collaborators at Cornell.  This version is based on release 1.0 of the Unix
   version. 

PCx is an interior-point predictor-corrector linear programming package and it
oficial home page is http://pages.cs.wisc.edu/~swright/PCx/. The
code has been developed at the Optimization Technology Center, which is a
collaboration between Argonne National Laboratory and Northwestern University.

PCx is designed as a stand-alone solver. Because of its modular structure and
fairly transparent data structures, it is not too difficult to integrate into
your application. Together with some of our users, we have recently investigated
new features such as a callable library and a MATLAB interface.

The official documentation of PCx is a Technical Report write by the authors
of PCx (:download:`dvi version <ref/PCx-user.dvi>` and :download:`ps version
<ref/PCx-user.ps>`).  This is a extension of the official documentation.

The source code and documentation for PCx can be obtained through the World Wide
Web in https://github.com/r-gaia-cs/PCx. If you need help try `this mail list
<https://www.listas.unicamp.br/mailman/listinfo/pcx-dev-l>`_.

Users:

.. toctree::
   :maxdepth: 2

   authors.rst
   license.rst
   install.rst
   invoking.rst
   specifications_file.rst
   samples.rst
   benchmark.rst

Developers:

.. toctree::
   :maxdepth: 2

   formulation.rst
   algorithm.rst
   linear_algebra.rst
   presolver.rst
   data-structures.rst
   data-structures2.rst
   debugging.rst
   reference.rst

C API:

.. toctree::
   :maxdepth: 2

   c_api/algorithm.rst
   c_api/basics.rst
   c_api/cblas.rst
   c_api/dcolumns.rst
   c_api/hash.rst
   c_api/io.rst
   c_api/jair.rst
   c_api/lpmps.rst
   c_api/main.rst
   c_api/memory.rst
   c_api/Ng-Peyton.rst
   c_api/parameters.rst
   c_api/PCx.rst
   c_api/pre.rst
   c_api/presolve.rst
   c_api/rcm.rst
   c_api/readmps.rst
   c_api/scale.rst
   c_api/solver.rst
   c_api/solve.rst
   c_api/split.rst
   c_api/timers.rst
   c_api/wrappers.rst
   c_api/writemps.rst
   c_api/wssmp.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

