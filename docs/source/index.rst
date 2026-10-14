.. toctree::
   :hidden:
   :caption: stko
   :maxdepth: 1

   Video Tutorials <video_tutorials>
   Molecular <molecular>
   Calculators <calculators>
   Optimizers <optimizers>
   Cage optimisation workflow <cage_optimisation>
   Host-guest optimisation <hg_optimisation>
   Using the topology for GFN-FF optimisation <using_topology_for_opt>
   Cage analysis <cage_analysis>
   Helpers <helpers>


.. toctree::
   :hidden:
   :caption: Modules
   :maxdepth: 1

   Modules <modules.rst>

Welcome to stko's documentation!
================================

| GitHub: https://github.com/JelfsMaterialsGroup/stko
| Discord: https://discord.gg/zbCUzuxe2B

.. tip::

  ⭐ Star us on `GitHub <https://github.com/JelfsMaterialsGroup/stko>`_! ⭐

.. figure:: _static/logo.png

Overview
========

`stko <https://github.com/JelfsMaterialsGroup/stko>`_ is a Python library for
performing optimizations and calculations on complex molecules built using
`stk <https://github.com/lukasturcani/stk>`_. In the case of
optimizations, a clone of :class:`stk.Molecule` is returned. For
calculators, a ``Results`` class are used to calculate and extract
properties of an :class:`stk.Molecule`.

Installation
============

:mod:`.stko` can be installed directly with pip:

.. code-block:: bash

  pip install stko

Some optional dependencies are only available through conda:

.. code-block:: bash

  # for gfn-ff with topology constrains
  # This will not work on Windows at this stage
  # (https://pprcht.github.io/gfnff/python.html)
  pip install gfnff

  # for xtb
  mamba install xtb
  # for openbabel, assuming you are not using Python >= 3.13!
  mamba install openbabel

With OpenMM
-----------

To get :mod:`.stko` and use ``OpenMM``, we had some installation issues. The
current solution is to first, in a new environment, install the ``OpenMM``
requirements:

.. code-block:: bash

  mamba install -c conda-forge openff-toolkit

Then install :mod:`.stko` with pip, but with the cuda variant to take advantage
of GPU speed up (note that this is a heavy installation!).

.. code-block:: bash

  pip install stko[cuda]

We also removed the default installation of ``espaloma_charge`` that provides
the ML-based ``espaloma-am1bcc`` partial charges method. If users want this
package, create a new environment and install their dependancies (if this
fails, please check their
`instructions <https://github.com/choderalab/espaloma-charge>`_), then install
``stko``:

.. code-block:: bash

    mamba install -c conda-forge espaloma_charge openff-toolkit
    pip install stko[cuda]


Developer Setup
---------------

#. Install `just`_.
#. In a new virtual environment run::

    $ just dev

#. Run code checks::

    $ just check

.. _`just`: https://github.com/casey/just

Dependencies
------------

The software packages we offer optimizers for are also depencies depending
on the desired functions used. These are:

* `MacroModel <https://www.schrodinger.com/platform/products/macromodel/>`_
* `GULP <http://gulp.curtin.edu.au/gulp/>`_
* `XTB <https://xtb-docs.readthedocs.io/en/latest/>`_
* `OpenBabel <https://github.com/openbabel/openbabel>`_
* `OpenMM <https://openmm.org/>`_
* `OpenFF <https://openforcefield.org/>`_


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

How To Contribute
-----------------

If you have any questions or find problems with the code, please submit
an issue.

If you wish to add your own code to this repository, please send us a
Pull Request. Please maintain the testing and style that is used
throughout ```stko``.


How To Cite
-----------

If you use ``stko`` please cite

    https://github.com/JelfsMaterialsGroup/stko
