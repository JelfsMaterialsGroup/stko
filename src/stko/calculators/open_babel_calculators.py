"""
OpenBabel Calculators
=====================

#. :class:`.OpenBabelEnergy`

Wrappers for calculators within the `openbabel` code.

"""

import logging
import os

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None

from stko.calculators.calculators import Calculator
from stko.calculators.results.energy_results import EnergyResults
from stko.utilities.exceptions import (
    ForceFieldSetupError,
    WrapperNotInstalledError,
)

logger = logging.getLogger(__name__)


class OpenBabelEnergy(Calculator):
    """
    Uses OpenBabel to calculate forcefield energies.[1]_

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose energy we want to know.
        mol1 = stk.BuildingBlock('CCCNCCCN')

        # Create the energy calculator.
        energy_calc = stko.OpenBabelEnergy('uff')

        # Calculate the energy.
        results = energy_calc.get_results(mol1)
        energy = results.get_energy()
        unit_string = results.get_unit_string()

    References
    ----------
    .. [1] http://openbabel.org/dev-api/classOpenBabel_1_
    1OBForceField.shtml#a2f2732698efde5c2f155bfac08fd9ded

    """

    def __init__(self, forcefield):
        """
        Initialize `openbabel` forcefield energy calculation.

        Parameters
        ----------
        forcefield : :class:`str`
            Forcefield to use. Options include `uff`, `gaff`,
            `ghemical`, `mmff94`.

        """

        if openbabel is None:
            raise WrapperNotInstalledError(
                "openbabel is not installed; see README for " "installation."
            )

        self._forcefield = forcefield

    def calculate(self, mol):
        temp_file = "temp.mol"
        mol.write(temp_file)
        try:
            obConversion = openbabel.OBConversion()
            obConversion.SetInFormat("mol")
            OBMol = openbabel.OBMol()
            obConversion.ReadFile(OBMol, temp_file)
            OBMol.PerceiveBondOrders()
        finally:
            os.system("rm temp.mol")

        forcefield = openbabel.OBForceField.FindForceField(self._forcefield)
        outcome = forcefield.Setup(OBMol)
        if not outcome:
            raise ForceFieldSetupError(
                f"Openbabel: {self._forcefield} could not be setup for {mol}"
            )

        yield forcefield.Energy()

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
            unit_string="kJ mol-1",
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
