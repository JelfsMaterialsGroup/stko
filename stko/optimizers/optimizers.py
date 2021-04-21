"""
Optimizers
==========

#. :class:`.NullOptimizer`
#. :class:`.OptimizerSequence`
#. :class:`.TryCatchOptimizer`
#. :class:`.RaisingOptimizer`


Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and can optimize a molecule
with :meth:`~.Optimizer.optimize`.

.. code-block:: python

    import stk
    import stko

    mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
    mmff = stko.MMFF()
    mol = mmff.optimize(mol)

    # Optimizers also work with ConstructedMolecule objects.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(mol, ),
            repeating_unit='A',
            num_repeating_units=3,
        )
    )
    etkdg = stko.ETKDG()
    polymer = etkdg.optimize(polymer)

Sometimes it is desirable to chain multiple optimizations, one after
another. For example, before running an optimization, it may be
desirable to embed a molecule first, to generate an initial structure.
:class:`.OptimizerSequence` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    optimizer_sequence = stko.OptimizerSequence(etkdg, mmff)

    # Run each optimizer in sequence.
    polymer = optimizer_sequence.optimize(polymer)

.. _`adding optimizers`:

Making New Optimizers
---------------------

New optimizers can be made by simply making a class which inherits the
:class:`.Optimizer` class. In addition to this, the new class must
define a :meth:`~.Optimizer.optimize` method. The method must take 1
mandatory `mol` parameter. :meth:`~.Optimizer.optimize` will take the
`mol` and change its structure in whatever way it likes. Beyond this
there are no requirements. New optimizers can be added into the
:mod:`.optimizers` submodule or into a new submodule.

"""

import logging


logger = logging.getLogger(__name__)


class Optimizer:
    """
    A base class for optimizers.

    """

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        """

        ...


class OptimizerSequence(Optimizer):
    """
    Applies optimizers in sequence.

    Examples
    --------
    Let's say we want to embed a molecule with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
        optimizer = stko.OptimizerSequence(stko.ETKDG(), stko.MMFF())
        mol = optimizer.optimize(mol)

    """

    def __init__(self, *optimizers):
        """
        Initialize a :class:`OptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

        """

        self._optimizers = optimizers

    def optimize(self, mol):

        for optimizer in self._optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol}".')
            mol = optimizer.optimize(mol)

        return mol


class TryCatchOptimizer(Optimizer):
    """
    Try to optimize with a :class:`Optimizer`, use another on failure.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        # Create some molecules to optimize.
        mol1 = stk.BuildingBlock('NCCN')
        mol2 = stk.BuildingBlock('CCCCC')
        mol3 = stk.BuildingBlock('O=CCCN')

        # Create an optimizer which may fail.
        uff = stko.UFF()

        # Create a backup optimizer.
        mmff = stko.MMFF()

        # Make an optimizer which tries to run raiser and if that
        # raises an error, will run mmff on the molecule instead.
        try_catch = stko.TryCatchOptimizer(
            try_optimizer=uff,
            catch_optimizer=mmff,
        )

        # Optimize the molecules. In each case if the optimization with
        # UFF fails, MMFF is used to optimize the molecule instead.
        mol1 = try_catch.optimize(mol1)
        mol2 = try_catch.optimize(mol2)
        mol3 = try_catch.optimize(mol3)

    """

    def __init__(self, try_optimizer, catch_optimizer):
        """
        Initialize a :class:`TryCatchOptimizer` instance.

        Parameters
        ----------
        try_optimizer : :class:`Optimizer`
            The optimizer which is used initially to try and optimize a
            :class:`.Molecule`.

        catch_optimizer : :class:`Optimizer`
            If `try_optimizer` raises an error, this optimizer is
            run on the :class:`.Molecule` instead.

        """

        self._try_optimizer = try_optimizer
        self._catch_optimizer = catch_optimizer

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        """

        try:
            return self._try_optimizer.optimize(mol)
        except Exception:
            try_name = self._try_optimizer.__class__.__name__
            catch_name = self._catch_optimizer.__class__.__name__
            logger.error(
                f'{try_name} failed, trying {catch_name}.',
                exc_info=True
            )
            return self._catch_optimizer.optimize(mol)
