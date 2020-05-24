"""
RDKit Calculators
=================

#. :class:`.MMFFEnergy`
#. :class:`.UFFEnergy`

Wrappers for calculators within the :mod:`rdkit` code.

"""

import logging
from rdkit.Chem import AllChem as rdkit

from ..base_calculator import Calculator


logger = logging.getLogger(__name__)


class MMFFEnergy(Calculator):
    """
    Uses the MMFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        mmff = stk.MMFFEnergy()

        # Calculate the energy.
        energy1 = mmff.get_energy(mol1)

    """

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.MMFFGetMoleculeForceField(
            rdkit_mol,
            rdkit.MMFFGetMoleculeProperties(rdkit_mol)
        )
        return ff.CalcEnergy()


class UFFEnergy(Calculator):
    """
    Uses the UFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        uff = stk.UFFEnergy()

        # Calculate the energy.
        energy1 = uff.get_energy(mol1)

    """

    def get_energy(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`float`
            The energy.

        """

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        # RingInfo needs to be initialized, else rdkit may raise an
        # error.
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.UFFGetMoleculeForceField(rdkit_mol)
        return ff.CalcEnergy()
