"""
Orca Calculators
================

#. :class:`.OrcaEnergy`

Wrappers for calculators within the :mod:`orca` code.

"""

import logging
import os
import shutil
import glob
import uuid
import subprocess as sp

from .calculators import Calculator
from .results import OrcaResults

logger = logging.getLogger(__name__)


class OrcaOptimizerError(Exception):
    ...


class OrcaEnergy(Calculator):
    """
    Uses Orca [1]_ to calculate energy and other properties.

    By default, :meth:`get_results` will extract other properties of
    the :class:`.Molecule` passed to :meth:`calculate`, which
    will be saved in the attributes of :class:`.OrcaResults`.

    All intermediate and output files from Orca are deleted at the end
    of the job (i.e. the `.gbw` file will be deleted) because they can
    quickly build up to large sizes. The `discard_output` option allows
    you to keep output files if desired.Additionally, the
    `write_input_only` option is available for jobs where you would
    like more customization or to run outside of the Python
    environment.

    Notes
    -----
    When running :meth:`calculate`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.OrcaEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Examples
    --------

    .. code-block:: python

        import stk
        import stko

        bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
        polymer = stk.ConstructedMolecule(
            stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit="AB",
                orientations=[0, 0],
                num_repeating_units=1
            )
        )

        # Optimize the constructed molecule so that it has a
        # reasonable structure.
        opt = stko.UFF()
        polymer = opt.optimize(polymer)

        # Calculate energy using Orca.
        orca = stko.OrcaEnergy(
            orca_path='/opt/orca/orca',
            topline='! SP B97-3c',
        )

        orca_results = orca.get_results(polymer)

        # Extract properties from the energy calculator for a given
        # molecule.
        total_energy = orca_results.get_total_energy()

    If you want the input file written (instead of the job run), you
    can use the `write_input_only` argument to save the input file
    in the `output_dir` as `orca_input.inp` with the input xyz file as
    `input_structure.xyz`.

    .. code-block:: python

        # Optimize the constructed molecule so that it has a
        # reasonable structure.
        optimizer = stko.ETKDG()
        polymer = optimizer.optimize(polymer)

        # Calculate energy using Orca.
        orca = stko.OrcaEnergy(
            orca_path='/opt/orca/orca',
            topline='! SP B97-3c',
            write_input_only=True,
        )

        orca.get_results(polymer)

    References
    ----------
    .. [1] https://orcaforum.kofo.mpg.de/app.php/portal

    """
    def __init__(
        self,
        orca_path,
        topline,
        basename=None,
        output_dir=None,
        num_cores=1,
        charge=0,
        multiplicity=1,
        write_input_only=False,
        discard_output=True,
    ):
        """
        Initializes a :class:`OrcaEnergy` instance.

        Parameters
        ----------
        orca_path : :class:`str`
            The path to the Orca executable.

        topline : :class:`str`
            Top line designating the type of calculation. Should start
            with `! `.

        basename : :class:`str`, optional
            Base name of Orca output files.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`, optional
            The number of cores Orca should use.

        charge : :class:`int`, optional
            Formal molecular charge.

        multiplicity : :class:`int`, optional
            Multiplicity of system (2S+1), where S is the spin.

        write_input_only : :class:`bool`, optional
            `True` if you only want the input file written and to not
            have the Orca job run.

        discard_output : :class:`bool`, optional
            `True` if you want to delete auxillary Orca output files
            such as the `.gbw` file.

        """

        self._orca_path = orca_path
        if basename is None:
            self._basename = f'_{str(uuid.uuid4().int)}'
        else:
            self._basename = basename
        self._output_dir = output_dir
        self._topline = topline
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._multiplicity = multiplicity
        self._write_input_only = write_input_only
        self._discard_output = discard_output

    def _write_input_file(self, path, xyz_file):
        # Write top line and base name.
        string = f'{self._topline}\n\n%base "{self._basename}"\n'

        # Add multiprocessing section.
        string += (
            f'%pal\n   nprocs {self._num_cores}\nend\n\n'
            f'%scf\n   MaxIter 2000\nend\n\n'
        )
        # Add geometry section.
        string += (
            f'* xyzfile {self._charge} {self._multiplicity} '
            f'{xyz_file}\n'
        )

        with open(path, 'w') as f:
            f.write(string)

    def _check_outcome(self, out_file):

        if not os.path.exists(out_file):
            raise OrcaOptimizerError(
                f'{out_file} does not exist, suggesting the job did '
                'not run.'
            )

        with open(out_file, 'r') as f:
            lines = f.readlines()
            if '****ORCA TERMINATED NORMALLY****' not in lines[-2]:
                raise OrcaOptimizerError(
                    'Orca job did not terminate normally.'
                )

        tmp_files = glob.glob(f'{self._basename}*tmp')
        if len(tmp_files) > 0:
            raise OrcaOptimizerError(
                'tmp files exist, suggesting the job did not complete '
                'or did not converge.'
            )

    def _clean_up(self):
        for to_del in glob.glob(f'{self._basename}*'):
            os.remove(to_del)

    def _run_orca(
        self,
        xyz_file,
        input_file,
        out_file,
        init_dir,
        output_dir,
    ):
        """
        Runs Orca.

        Parameters
        ----------
        xyz_file : :class:`str`
            The name of the input structure ``.xyz`` file.

        input_file : :class:`str`
            The name of input file to be written.

        out_file : :class:`str`
            The name of output file with Orca results.

        init_dir : :class:`str`
            The name of the current working directory.

        output_dir : :class:`str`
            The name of the directory into which files generated during
            the calculation are written.

        Returns
        -------
        None : :class:`NoneType`

        """

        cmd = f'{self._orca_path} {input_file}'

        try:
            os.chdir(output_dir)
            self._write_input_file(input_file, xyz_file)
            if not self._write_input_only:
                with open(out_file, 'w') as f:
                    # Note that sp.call will hold the program until
                    # completion of the calculation.
                    sp.call(
                        cmd,
                        stdin=sp.PIPE,
                        stdout=f,
                        stderr=sp.PIPE,
                        # Shell is required to run complex arguments.
                        shell=True
                    )
                self._check_outcome(out_file)
                if self._discard_output:
                    self._clean_up()
        finally:
            os.chdir(init_dir)

    def calculate(self, mol):
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)

        init_dir = os.getcwd()
        xyz_file = os.path.join(output_dir, 'input_structure.xyz')
        input_file = os.path.join(output_dir, 'orca_input.inp')
        out_file = os.path.join(output_dir, 'orca_energy.output')
        mol.write(xyz_file)
        yield self._run_orca(
            xyz_file=xyz_file,
            input_file=input_file,
            out_file=out_file,
            init_dir=init_dir,
            output_dir=output_dir,
        )

    def get_results(self, mol):
        """
        Calculate the Orca properties of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy is to be calculated.

        Returns
        -------
        :class:`.OrcaResults`
            The properties, with units, from Orca calculations.

        """

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        out_file = os.path.join(output_dir, 'orca_energy.output')

        if self._write_input_only:
            next(self.calculate(mol))
        else:
            return OrcaResults(
                generator=self.calculate(mol),
                output_file=out_file,
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

        return self.get_results(mol).get_total_energy()[0]
