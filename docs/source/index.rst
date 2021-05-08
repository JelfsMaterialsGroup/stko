Welcome to stko's documentation!
================================

============
Introduction
============

GitHub: https://github.com/JelfsMaterialsGroup/stko



Overview
========

`stko <https://github.com/lukasturcani/stk>`_ is a Python library which performs optimizations and calculations on
complex molecules built using `stk <https://github.com/JelfsMaterialsGroup/stko>`_.
Optimizations or calculations are typically performed by providing an :class:`stk.ConstructedMolecule` to the function.
In the case of optimizations, a clone of :class:`stk.ConstructedMolecule` is returned
and calculators return a :class:`float`.


Documentation
=============

.. glossary::

   Version
      |version|


.. toctree::
   :hidden:
   :caption: Installation
   :maxdepth: 1

   Install <install>


.. toctree::
   :hidden:
   :caption: Optimizers
   :maxdepth: 2

   XTB <stko.optimizers.xtb>
   RDKit <stko.optimizers.rdkit>
   GULP <stko.optimizers.gulp>
   Collapser <stko.optimizers.collapser>
   MacroModel <stko.optimizers.macromodel>


.. toctree::
   :hidden:
   :caption: Calculators
   :maxdepth: 2

   RDKit <stko.calculators.rdkit_calculators>
   XTB <stko.calculators.xtb_calculators>


.. toctree::
   :hidden:
   :caption: Molecular Systems
   :maxdepth: 2

   Topology Extraction <stko.molecular.topology_extractor>
   Unit Cell <stko.molecular.periodic.unitcell>
   Network Material <stko.molecular.networkx>


Indices and Tables
==================

* :ref:`genindex`
* :ref:`py-modindex`