"""
RDKit Optimizers
================

#. :class:`.MMFF`
#. :class:`.UFF`
#. :class:`.ETKDG`

Wrappers for optimizers within the :mod:`rdkit` code.

.. code-block:: python

    import stk
    import stko

    # Optimizers work on stk.Molecule.
    polymer = stk.ConstructedMolecule(
        building_blocks=[mol],
        topology_graph=stk.polymer.Linear('A', [0], n=3)
    )
    etkdg = stk.ETKDG()
    etkdg.optimize(polymer)

"""

import logging
import rdkit.Chem.AllChem as rdkit

from .optimizer import Optimizer


logger = logging.getLogger(__name__)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        mmff = stk.MMFF()
        mmff.optimize(mol)

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
        None : :class:`NoneType`

        """

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        uff = stk.UFF()
        uff.optimize(mol)

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
        None : :class:`NoneType`

        """

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        etkdg = stk.ETKDG()
        etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12, use_cache=False):
        """
        Initialize a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._random_seed = random_seed
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        mol.set_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )
