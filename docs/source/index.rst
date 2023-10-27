.. toctree::
   :hidden:
   :caption: stko
   :maxdepth: 2

   Molecular <molecular>
   Calculators <calculators>
   Optimizers <optimizers>


.. toctree::
   :hidden:
   :caption: Modules
   :maxdepth: 1

   Modules <modules.rst>

Welcome to stko's documentation!
================================

GitHub: https://github.com/JelfsMaterialsGroup/stko


Install
=======

:mod:`.stko` can be installed directly with pip::

   $ pip install stko


Dependencies
------------

The software packages we offer optimizers are also depencies depending
on the desired functions used. These are:

* `MacroModel <https://sites.google.com/site/orcainputlibrary/home/>`_
* `GULP <http://gulp.curtin.edu.au/gulp/>`_
* `XTB <https://xtb-docs.readthedocs.io/en/latest/>`_
* `OpenBabel <https://github.com/openbabel/openbabel>`_


Overview
========

`stko <https://github.com/lukasturcani/stk>`_ is a Python library which
performs optimizations and calculations on complex molecules built using
`stk <https://github.com/JelfsMaterialsGroup/stko>`_. In the case of
optimizations, a clone of :class:`stk.Molecule` is returned. For
calculators, a :class:`.Calculator` and a :class:`.Results` are used to
calculate and extract properties of an :class:`stk.Molecule`.


Indices and Tables
==================

* :ref:`genindex`
