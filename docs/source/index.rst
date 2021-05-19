Welcome to stko's documentation!
================================

GitHub: https://github.com/JelfsMaterialsGroup/stko


Overview
========

`stko <https://github.com/lukasturcani/stk>`_ is a Python library which performs optimizations and calculations on
complex molecules built using `stk <https://github.com/JelfsMaterialsGroup/stko>`_.
Optimizations or calculations are typically performed by providing an :class:`stk.ConstructedMolecule` to the function.
In the case of optimizations, a clone of :class:`stk.ConstructedMolecule` is returned.
For calculators, a :class:`.Calculator` and a :class:`.Results` are used to calculate and extract properties of an :class:`stk.Molecule`.


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

   Optimizers <stko.optimizers.optimizers>
   XTB <stko.optimizers.xtb>
   RDKit <stko.optimizers.rdkit>
   GULP <stko.optimizers.gulp>
   Collapser <stko.optimizers.collapser>
   MacroModel <stko.optimizers.macromodel>
   
   
.. toctree::
   :hidden:
   :caption: Calculators
   :maxdepth: 2

   Calculators <stko.calculators.calculators>
   RDKit <stko.calculators.rdkit_calculators>
   XTB <stko.calculators.xtb_calculators>
   Shape <stko.calculators.shape_calculators>
   Torsion <stko.calculators.torsion_calculators>


.. toctree::
   :hidden:
   :caption: Results
   :maxdepth: 2

   Results <stko.calculators.results.results>
   Energy Results <stko.calculators.results.energy_results>
   XTB Results <stko.calculators.results.xtb_results>
   Shape Results <stko.calculators.results.shape_results>
   Torsion Results <stko.calculators.results.torsion_results>


.. toctree::
   :hidden:
   :caption: Molecular Systems
   :maxdepth: 2

   Topology Extraction <stko.molecular.topology_extractor>
   Unit Cell <stko.molecular.periodic.unitcell>
   NetworkX <stko.molecular.networkx>
   Torsions <stko.molecular.torsion>


Indices and Tables
==================

* :ref:`genindex`
