"""
MacorModel Calculators
=================

# . :class:`.MacroModelForceFieldEnergy`

Wrappers for calculators that use Schrödinger's MacroModel software.

"""

import logging
from .calculators import Calculator
from ..optimizers import MacroModel
from .results import EnergyResults
from uuid import uuid4
from ..utilities import move_generated_macromodel_files
import re

logger = logging.getLogger(__name__)


class MacroModelCalculator(MacroModel, Calculator):
    """
    Base class for MacroModel calculators

    """

    def __init__(
        self,
        macromodel_path,
        output_dir,
        timeout,
        force_field,
    ):
        MacroModel.__init__(
            self,
            macromodel_path=macromodel_path,
            output_dir=output_dir,
            timeout=timeout,
            force_field=force_field,
            # Below argument are unnecesary for calculators.
            # Needed for MacroModel.__init__
            maximum_iterations=2500,
            minimum_gradient=0.05
        )
        Calculator.__init__(self)


class MacroModellForceFieldEnergy(MacroModelCalculator):
    """
    Uses MacroModel software to calculate energies.

    Examples
    --------
    .. code-block:: python

        import stk
        import stko

        # Create a molecule whose energy
        we want to know
        mol = stk.BuildingBlock('NCCN')

        # Create energy calculator
        opls_ff = stko.MacroModellForceFieldEnergy(#Insert MacroModel
        path here)

        # Calculate the energy
        results = opls_ff.get_results(mol)
        energy = opls_ff.get_energy(mol)
        unit_string = results.get_unit_string()

    """

    def __init__(self,
                 macromodel_path,
                 output_dir=None,
                 timeout=None,
                 force_field=16,
                 ):
        """
        Initilize a :class:`.MacroModellForceFieldEnergy` instance.

        Parameters
        ----------
        macromodel_path : :class:`str`
            The full path to the MacroModel executable.

        output_dir : :class:`str`
            The path to the output directory.

        timeout : :class:`int`
            The number of seconds to wait for the calculation to
            complete.

        forcefield : :class:`int`
            The force field to use.
            Force field arguments can be the following:
            +------------+------------+
            |  1  | MM2               |
            +------------+------------+
            |  2  | MM3               |
            +------------+------------+
            |  3  | AMBER             |
            +------------+------------+
            |  4  | AMBER94           |
            +------------+------------+
            |  5  | OPLSA             |
            +------------+------------+
            |  10 | MMFF94 and MMFF94s|
            +------------+------------+
            |  14 | OPLS_2005         |
            +------------+------------+
            |  16 | OPLS3e            |
            +------------+------------+

        """
        super().__init__(
            macromodel_path=macromodel_path, output_dir=output_dir,
            timeout=timeout,
            force_field=force_field,
        )
        self._macromodel_path = macromodel_path
        self._output_dir = output_dir
        self._timeout = timeout
        self._forcefield = force_field

    def get_results(self, mol):
        """
        Calculate the energy of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to calculate the energy of.

        Returns
        -------
        :class:`.EnergyResults`
            The energy and units of the energy.

        """

        return EnergyResults(
            generator=self.calculate(mol),
            unit_string='kJ mol-1'
        )

    def calculate(self, mol):
        """
        Performs the calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to calculate the energy of.

        """
        run_name = str(uuid4().int)
        if self._output_dir is None:
            output_dir = run_name
        else:
            output_dir = self._output_dir
        mol_path = f'{run_name}.mol'
        mae_path = f'{run_name}.mae'
        # First write the molecule
        mol.write(mol_path)
        # MacroModel requires a `.mae` file
        self._run_structconvert(mol_path, mae_path)
        # Generate the `.com` file
        self._generate_com(mol, run_name)
        # Run the calculation
        self._run_bmin(mol, run_name)
        # Convert `.maegz` output to `.mae`
        self._convert_maegz_to_mae(run_name)
        # Read the results
        result = self._read_mmo(run_name)
        move_generated_macromodel_files(run_name, output_dir)
        yield result

    @staticmethod
    def _read_mmo(run_name):
        """
        Reads and parses the `.mmo` file following the calculation.

        Parameters
        ----------
        run_name : :class:`str`
            The name of the calculation.

        Returns
        -------
        :class:`.MacroModelOutput`
            The results of the calculation.

        """

        mmo_path = f'{run_name}-out.mmo'
        with open(mmo_path, 'r') as mmo:
            mmo_lines = mmo.readlines()
        # Read the energy
        e_line = list(
            filter(
                lambda line: 'Total energy' in
                line,
                mmo_lines
            )
        )

        return float(re.findall(r'[-+]?\d*\.\d+|\d+', e_line[0])[0])

    def _generate_com(self, mol, run_name):
        """
        Create a `.com` file for a MacroModel energy calculation.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule to calculate the energy of.

        run_name: : class: `str`
            The name of the calculation. Files
            generated by MacroModel will have this name.

        Returns
        -------
        None: : class: `NoneType`

        """

        logger.debug(f'Generating .com file for "{mol}".')
        # Create the body of the `.com` file
        # Select the force field
        line1 = ('FFLD', self._force_field, 1, 0, 0, 1, 0, 0, 0)
        line2 = ('READ', 0, 0, 0, 0, 0, 0, 0, 0)
        # Line 3 writes the complete energy to the output log file
        line3 = ('ELST', 0, 0, 0, 0, 0, 0, 0, 0)

        com_block = '\n'.join([
            self._get_com_line(*line1),
            self._get_com_line(*line2),
            self._get_com_line(*line3),
        ])

        with open(f'{run_name}.com', 'w') as com:
            # Name of the MacroModel file
            com.write(f'{run_name}.mae\n')
            # Name of the output file
            com.write(f'{run_name}-out.maegz\n')
            com.write(com_block)
