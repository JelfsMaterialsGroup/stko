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

    mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])

    # Optimizers work on stk.Molecule.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(mol, ),
            repeating_unit='A',
            num_repeating_units=3,
        )
    )
    etkdg = stko.ETKDG()
    polymer = etkdg.optimize(polymer)

"""

import logging
import rdkit.Chem.AllChem as rdkit

from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCNCCN')
        mmff = stko.MMFF()
        mol = mmff.optimize(mol)

    """

    def __init__(self, ignore_inter_interactions=True):

        self._ignore_inter_interactions = (
            ignore_inter_interactions
        )

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(
            rdkit_mol,
            ignoreInterfragInteractions=self._ignore_inter_interactions
        )
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )

        return mol


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCNCCN')
        uff = stko.UFF()
        mol = uff.optimize(mol)

    """

    def __init__(self, ignore_inter_interactions=True):

        self._ignore_inter_interactions = ignore_inter_interactions

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(
            rdkit_mol,
            ignoreInterfragInteractions=self._ignore_inter_interactions
        )
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )

        return mol


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCNCCN')
        etkdg = stko.ETKDG()
        mol = etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12):
        """
        Initialize a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._random_seed = random_seed

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

        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )

        return mol
