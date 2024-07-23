import logging
import os
import shutil
import subprocess as sp
import uuid
from collections import abc
from pathlib import Path

import stk

from stko._internal.calculators.results.orca_results import OrcaResults
from stko._internal.utilities.exceptions import OptimizerError, PathError

logger = logging.getLogger(__name__)


class OrcaEnergy:
    """Uses Orca to calculate energy and other properties.

    By default, :meth:`get_results` will extract other properties of
    the :class:`stk.Molecule` passed to :meth:`calculate`, which
    will be saved in the attributes of :class:`stko.OrcaResults`.

    All intermediate and output files from Orca are deleted at the end
    of the job (i.e. the ``.gbw`` file will be deleted) because they can
    quickly build up to large sizes. The `discard_output` option allows
    you to keep output files if desired.Additionally, the
    `write_input_only` option is available for jobs where you would
    like more customization or to run outside of the Python
    environment.

    See Also:
        * Orca: https://orcaforum.kofo.mpg.de/app.php/portal

    Parameters:
        orca_path:
            The path to the Orca executable.

        topline:
            Top line designating the type of calculation. Should start
            with ``!``.

        basename:
            Base name of Orca output files.

        output_dir:
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores:
            The number of cores Orca should use.

        charge:
            Formal molecular charge.

        multiplicity:
            Multiplicity of system (2S+1), where S is the spin.

        write_input_only:
            ``True`` if you only want the input file written and to not
            have the Orca job run.

        discard_output:
            ``True`` if you want to delete auxillary Orca output files
            such as the ``.gbw`` file.

    Notes:
        When running :meth:`calculate`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        If multi-threading is being used an error could occur if two
        different threads need to know about the current working directory
        as :class:`stko.OrcaEnergy` can change it from under them.

        Note that this does not have any impact on multi-processing,
        which should always be safe.

    Examples:
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


    """

    def __init__(  # noqa: PLR0913
        self,
        orca_path: Path | str,
        topline: str,
        basename: str | None = None,
        output_dir: Path | str | None = None,
        num_cores: int = 1,
        charge: int = 0,
        multiplicity: int = 1,
        write_input_only: bool = False,
        discard_output: bool = True,
    ) -> None:
        orca_path = Path(orca_path)
        self._check_path(orca_path)
        self._orca_path = orca_path
        if basename is None:
            self._basename = f"_{uuid.uuid4().int!s}"
        else:
            self._basename = basename
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._topline = topline
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._multiplicity = multiplicity
        self._write_input_only = write_input_only
        self._discard_output = discard_output

    def _check_path(self, path: Path) -> None:
        if not path.exists():
            msg = f"Orca not found at {path}"
            raise PathError(msg)

    def _write_input_file(self, path: Path, xyz_file: Path) -> None:
        # Write top line and base name.
        string = f'{self._topline}\n\n%base "{self._basename}"\n'

        # Add multiprocessing section.
        string += (
            f"%pal\n   nprocs {self._num_cores}\nend\n\n"
            f"%scf\n   MaxIter 2000\nend\n\n"
        )
        # Add geometry section.
        string += (
            f"* xyzfile {self._charge} {self._multiplicity} " f"{xyz_file}\n"
        )

        path.write_text(string)

    def _check_outcome(self, out_file: Path) -> None:
        if not out_file.exists():
            msg = (
                f"ORCA: {out_file} does not exist, suggesting the job did "
                "not run."
            )
            raise OptimizerError(msg)

        with out_file.open() as f:
            lines = f.readlines()
            if "****ORCA TERMINATED NORMALLY****" not in lines[-2]:
                msg = "ORCA: Orca job did not terminate normally."
                raise OptimizerError(msg)

        tmp_files = list(Path().glob(f"{self._basename}*tmp"))
        if len(tmp_files) > 0:
            msg = (
                "ORCA: tmp files exist, suggesting the job did not complete "
                "or did not converge."
            )
            raise OptimizerError(msg)

    def _clean_up(self) -> None:
        for to_del in Path().glob(f"{self._basename}*"):
            to_del.unlink()

    def _run_orca(
        self,
        xyz_file: Path,
        input_file: Path,
        out_file: Path,
        init_dir: Path,
        output_dir: Path,
    ) -> None:
        """Runs Orca.

        Parameters:
            xyz_file:
                The name of the input structure ``.xyz`` file.

            input_file:
                The name of input file to be written.

            out_file:
                The name of output file with Orca results.

            init_dir:
                The name of the current working directory.

            output_dir:
                The name of the directory into which files generated during
                the calculation are written.

        """
        cmd = f"{self._orca_path} {input_file}"

        try:
            os.chdir(output_dir)
            self._write_input_file(input_file, xyz_file)
            if not self._write_input_only:
                with out_file.open("w") as f:
                    # Note that sp.call will hold the program until
                    # completion of the calculation.
                    sp.call(  # noqa: S602
                        cmd,
                        stdin=sp.PIPE,
                        stdout=f,
                        stderr=sp.PIPE,
                        # Shell is required to run complex arguments.
                        shell=True,
                    )
                self._check_outcome(out_file)
                if self._discard_output:
                    self._clean_up()
        finally:
            os.chdir(init_dir)

    def calculate(self, mol: stk.Molecule) -> abc.Generator:
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        init_dir = Path.cwd()
        xyz_file = output_dir / "input_structure.xyz"
        input_file = output_dir / "orca_input.inp"
        out_file = output_dir / "orca_energy.output"
        mol.write(xyz_file)
        self._run_orca(
            xyz_file=xyz_file,
            input_file=input_file,
            out_file=out_file,
            init_dir=init_dir,
            output_dir=output_dir,
        )
        yield

    def get_results(self, mol: stk.Molecule) -> OrcaResults | None:
        """Calculate the Orca properties of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The properties, with units, from Orca calculations or ``None``
            if ``write_input_only`` mode.

        """
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        out_file = output_dir / "orca_energy.output"

        if self._write_input_only:
            next(self.calculate(mol))
            return None
        return OrcaResults(
            generator=self.calculate(mol),
            output_file=out_file,
        )

    def get_energy(self, mol: stk.Molecule) -> float | None:
        """Calculate the energy of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The energy or ``None`` if ``write_input_only`` mode.

        """
        results = self.get_results(mol)
        if results is None:
            return None
        return results.get_total_energy()[0]
