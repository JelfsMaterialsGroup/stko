.. toctree::
   :hidden:
   :caption: stko
   :maxdepth: 1

   Video Tutorials <video_tutorials>
   Molecular <molecular>
   Calculators <calculators>
   Optimizers <optimizers>
   Cage analysis <cage_analysis>
   Helpers <helpers>


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

:mod:`.stko` can be installed directly with pip:

.. code-block:: bash

  pip install stko

Some optional dependencies are only available through conda:

.. code-block:: bash

  # for OpenMM
  mamba install openff-toolkit openmm openmmtools
  # for xtb
  mamba install xtb
  # for openbabel
  mamba install openbabel
  # for mdanalysis
  mamba install mdanalysis


Dependencies
------------

The software packages we offer optimizers are also depencies depending
on the desired functions used. These are:

* `MacroModel <https://sites.google.com/site/orcainputlibrary/home/>`_
* `GULP <http://gulp.curtin.edu.au/gulp/>`_
* `XTB <https://xtb-docs.readthedocs.io/en/latest/>`_
* `OpenBabel <https://github.com/openbabel/openbabel>`_
* `OpenMM <https://openmm.org/>`_
* `OpenFF <https://openforcefield.org/>`_


Overview
========

`stko <https://github.com/JelfsMaterialsGroup/stko>`_ is a Python library which
performs optimizations and calculations on complex molecules built using
`stk <https://github.com/lukasturcani/stk>`_. In the case of
optimizations, a clone of :class:`stk.Molecule` is returned. For
calculators, a ``Results`` class are used to calculate and extract
properties of an :class:`stk.Molecule`.

Examples
--------

For every class (including ``Calculator``, ``Optimizer``), there are small
examples of usage on the associated docs page. We have a page dedicated to
analysing `cage structures <cage_analysis.html>`_. There are also some examples
for ``stko`` usage available `here <https://github.com/JelfsMaterialsGroup/stko/tree/master/examples>`_.
These cover:

* `Basic examples <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/basic_example.py>`_
* `Molecule alignment <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/aligner_example.py>`_
* `Using calculators <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/calculators_example.py>`_
* `Cage analysis <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/cage_analysis_example.py>`_
* `Using Gulp <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/gulp_test_example.py>`_
* `Splitting molecules <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/molecule_splitter_example.py>`_
* `Interfacing with MDAnalysis <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/mdanalysis_example.py>`_
* `Interfacing with OpenBabel <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/obabel_example.py>`_
* `Interfacing with Orca <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/orca_example.py>`_
* `Calculating molecular shape with RDKit <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/shape_example.py>`_
* `Extracting stk topology graphs from molecules <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/topology_extraction_example.py>`_
* `Analysing torsions <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/torsion_example.py>`_
* `Converting molecules to their Zmatrix <https://github.com/JelfsMaterialsGroup/stko/blob/master/examples/zmatrix_example.py>`_

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
