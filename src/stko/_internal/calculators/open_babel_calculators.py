import logging
from collections import abc
from pathlib import Path

import stk

try:
    from openbabel import openbabel
except ImportError:
    openbabel = None


from stko._internal.calculators.results.energy_results import EnergyResults
from stko._internal.utilities.exceptions import (
    ForceFieldSetupError,
    WrapperNotInstalledError,
)

logger = logging.getLogger(__name__)


class OpenBabelEnergy:
    """Uses OpenBabel to calculate forcefield energies.

    Parameters:
        forcefield:
            Forcefield to use. Options include `uff`, `gaff`,
            `ghemical`, `mmff94`.

    Raises:
        :class:`WrapperNotInstalledError`: if `openbabel` not installed.


    See Also:
        * OpenBabel: https://github.com/openbabel/openbabel

    Examples:
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


    """

    def __init__(self, forcefield: str) -> None:
        if openbabel is None:
            msg = "openbabel is not installed; see README for installation."
            raise WrapperNotInstalledError(msg)

        self._forcefield = forcefield

    def calculate(self, mol: stk.Molecule) -> abc.Iterable[float]:
        temp_file = "temp.mol"
        mol.write(temp_file)
        try:
            ob_conversion = openbabel.OBConversion()
            ob_conversion.SetInFormat("mol")
            ob_mol = openbabel.OBMol()
            ob_conversion.ReadFile(ob_mol, temp_file)
            ob_mol.PerceiveBondOrders()
        finally:
            Path("temp.mol").unlink()

        forcefield = openbabel.OBForceField.FindForceField(self._forcefield)
        outcome = forcefield.Setup(ob_mol)
        if not outcome:
            msg = f"Openbabel: {self._forcefield} could not be setup for {mol}"
            raise ForceFieldSetupError(msg)

        yield forcefield.Energy()

    def get_results(self, mol: stk.Molecule) -> EnergyResults:
        """Calculate the energy of `mol`.

        Parameters:
            mol
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The energy and units of the energy.

        """
        return EnergyResults(
            generator=self.calculate(mol),
            unit_string="kJ mol-1",
        )

    def get_energy(self, mol: stk.Molecule) -> float:
        """Calculate the energy of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The energy.

        """
        return self.get_results(mol).get_energy()
