"""
RDKit Calculators
=================

#. :class:`.MMFFEnergy`
#. :class:`.UFFEnergy`

Wrappers for calculators within the :mod:`rdkit` code.

"""

import logging
from rdkit.Chem import AllChem as rdkit

from .calculators import Calculator
from .results import EnergyResults


logger = logging.getLogger(__name__)


class MMFFEnergy(Calculator):
    """
    Uses the MMFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko


        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        mmff = stko.MMFFEnergy()

        # Calculate the energy.
        results = mmff.get_results(mol1)
        energy = results.get_energy()
        unit_string = results.get_unit_string()

    """

    def __init__(self, ignore_inter_interactions=True):

        self._ignore_inter_interactions = (
            ignore_inter_interactions
        )

    def calculate(self, mol):
        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.MMFFGetMoleculeForceField(
            rdkit_mol,
            rdkit.MMFFGetMoleculeProperties(rdkit_mol),
            ignoreInterfragInteractions=self._ignore_inter_interactions
        )
        yield ff.CalcEnergy()

    def get_results(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.EnergyResults`
            The energy and units of the energy.

        """

        return EnergyResults(
            generator=self.calculate(mol),
            unit_string='kcal mol-1',
        )

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

        return self.get_results(mol).get_energy()


class UFFEnergy(Calculator):
    """
    Uses the UFF force field to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko


        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        uff = stko.UFFEnergy()

        # Calculate the energy.
        results = uff.get_results(mol1)
        energy = results.get_energy()
        unit_string = results.get_unit_string()

    """

    def __init__(self, ignore_inter_interactions=True):

        self._ignore_inter_interactions = (
            ignore_inter_interactions
        )

    def calculate(self, mol):
        rdkit_mol = mol.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        # RingInfo needs to be initialized, else rdkit may raise an
        # error.
        rdkit.GetSSSR(rdkit_mol)
        ff = rdkit.UFFGetMoleculeForceField(
            mol=rdkit_mol,
            ignoreInterfragInteractions=self._ignore_inter_interactions
        )
        yield ff.CalcEnergy()

    def get_results(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.EnergyResults`
            The energy and units of the energy.

        """

        return EnergyResults(
            generator=self.calculate(mol),
            unit_string='kcal mol-1',
        )

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

        return self.get_results(mol).get_energy()
