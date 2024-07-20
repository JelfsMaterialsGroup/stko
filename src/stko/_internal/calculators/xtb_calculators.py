import logging
import os
import shutil
import subprocess as sp
import uuid
from collections import abc
from pathlib import Path

import stk

from stko._internal.calculators.results.xtb_results import XTBResults
from stko._internal.utilities.exceptions import InvalidSolventError, PathError
from stko._internal.utilities.utilities import is_valid_xtb_solvent

logger = logging.getLogger(__name__)


class XTBEnergy:
    """Uses GFN-xTB to calculate energy and other properties.

    By default, :meth:`get_results` will extract other properties of
    the :class:`stk.Molecule` passed to :meth:`calculate`, which
    will be saved in the attributes of :class:`stko.XTBResults`.

    See Also:
        * GFN-xTB https://xtb-docs.readthedocs.io/en/latest/setup.html

    Parameters:
        xtb_path:
            The path to the xTB executable.

        gfn_version:
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir:
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores:
            The number of cores xTB should use.

        calculate_free_energy:
            Whether to calculate the total free energy and vibrational
            frequencies. Setting this to ``True`` can drastically
            increase calculation time and memory requirements.

        calculate_ip_and_ea:
            Whether to calculate the vertical ionisation potential and
            vertical electron affinity. Equivalent to  `--vipea` flag
            in command-line interface.

        electronic_temperature:
            Electronic temperature in Kelvin.

        solvent_model:
            Solvent model to use out of older `gbsa` and newer `alpb`.
            `gbsa` is default for backwards compatability, but `alpb`
            is recommended.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent:
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent_grid:
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        charge:
            Formal molecular charge.

        num_unpaired_electrons:
            Number of unpaired electrons.

        unlimited_memory:
            If ``True`` :meth:`energy` will be run without constraints
            on the stack size. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        write_sasa_info:
            If ``True``, the detailed info input will request ``gbsa=True`` and
            output SASA information from xtb. Requires a solvent model to be
            used.

    Notes:
        When running :meth:`calculate`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        If multi-threading is being used an error could occur if two
        different threads need to know about the current working directory
        as :class:`stko.XTBEnergy` can change it from under them.

        Note that this does not have any impact on multi-processing,
        which should always be safe.

        We thank Andrew Tarzia and Alejandro Santana-Bonilla for their
        contributions to this code.

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
            opt = stko.OptimizerSequence(
                stko.UFF(),
                stko.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
            )
            polymer = opt.optimize(polymer)

            # Calculate energy using GFN-xTB.
            xtb = stko.XTBEnergy(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True
            )

            xtb_results = xtb.get_results(polymer)

            # Extract properties from the energy calculator for a given
            # molecule.
            total_energy = xtb_results.get_total_energy()
            homo_lumo_gap = xtb_results.get_homo_lumo_gap()
            fermi_levels = xtb_results.get_fermi_level()
            homo_lumo_orbitals = xtb_results.get_homo_lumo_orbitals()
            full_dipole_moments = xtb_results.get_full_dipole_moments()

        If `calculate_free_energy` is ``True``, xTB performs a
        numerical Hessian calculation and calculates the total free energy
        and vibrational frequencies of a molecule. It is recommended that a
        well optimized structure be used as input for these calculations

        .. code-block:: python

            # Optimize the constructed molecule so that it has a
            # reasonable structure.
            optimizer = stko.OptimizerSequence(
                stko.ETKDG(),
                stko.XTB(
                    xtb_path='/opt/gfnxtb/xtb',
                    unlimited_memory=True,
                    opt_level='verytight'
                )
            )
            polymer = optimizer.optimize(polymer)

            # Calculate energy using GFN-xTB.
            xtb = stko.XTBEnergy(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True,
                calculate_free_energy=True
            )

            xtb_results = xtb.get_results(polymer)

            # Extract properties from the energy calculator for a given
            # molecule.
            total_free_energy = xtb_results.get_total_free_energy()
            total_frequencies = xtb_results.get_frequencies()

        If `calculate_ip_and_ea` is ``True``, xTB performs multiple single-
        point-energy calculations to calculate the vertical electron
        affinity and ionisation potential. It is recommended that a
        well optimized structure be used as input for these calculations

        .. code-block:: python

            # Optimize the constructed molecule so that it has a
            # reasonable structure.
            optimizer = stko.OptimizerSequence(
                stko.ETKDG(),
                stko.XTB(
                    xtb_path='/opt/gfnxtb/xtb',
                    unlimited_memory=True,
                    opt_level='verytight'
                )
            )
            polymer = optimizer.optimize(polymer)

            # Calculate energy using GFN-xTB.
            xtb = stko.XTBEnergy(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True,
                calculate_ip_and_ea=True,
            )

            xtb_results = xtb.get_results(polymer)

            # Extract properties from the energy calculator for a given
            # molecule.
            ip = xtb_results.get_ionisation_potential()
            ea = xtb_results.get_electron_affinity()

    """

    def __init__(  # noqa: PLR0913
        self,
        xtb_path: Path | str,
        gfn_version: int = 2,
        output_dir: Path | str | None = None,
        num_cores: int = 1,
        calculate_free_energy: bool = False,
        calculate_ip_and_ea: bool = False,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        solvent_grid: str = "normal",
        charge: int = 0,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,
        write_sasa_info: bool = False,
    ) -> None:
        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                msg = "XTB: No solvent valid for version", f" {gfn_version!r}."
                raise InvalidSolventError(msg)
            if not is_valid_xtb_solvent(
                gfn_version=gfn_version,
                solvent_model=solvent_model,
                solvent=solvent,
            ):
                msg = (
                    f"XTB: Solvent {solvent!r} and model {solvent_model!r}",
                    f" is invalid for version {gfn_version!r}.",
                )
                raise InvalidSolventError(msg)

        xtb_path = Path(xtb_path)
        self._check_path(xtb_path)
        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._num_cores = str(num_cores)
        self._calculate_free_energy = calculate_free_energy
        self._calculate_ip_and_ea = calculate_ip_and_ea
        self._electronic_temperature = str(electronic_temperature)
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self._solvent = solvent
        self._solvent_model = solvent_model
        self._solvent_grid = solvent_grid
        self._write_sasa_info = write_sasa_info

    def _check_path(self, path: Path) -> None:
        if not path.exists():
            msg = f"XTB not found at {path}"
            raise PathError(msg)

    def _write_detailed_control(self) -> None:
        sasa_info = "$write\n gbsa=true\n" if self._write_sasa_info else ""
        string = f"$gbsa\n gbsagrid={self._solvent_grid}\n{sasa_info}"

        Path("det_control.in").write_text(string)

    def _run_xtb(
        self,
        xyz: Path,
        out_file: Path,
        init_dir: Path,
        output_dir: Path,
    ) -> None:
        """Runs GFN-xTB.

        Parameters:
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

            init_dir:
                The name of the current working directory.

            output_dir:
                The name of the directory into which files generated during
                the calculation are written.

        """
        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        if self._solvent is not None:
            solvent = f"--{self._solvent_model} {self._solvent} "
        else:
            solvent = ""

        hess = "--hess" if self._calculate_free_energy else ""

        vipea = "--vipea" if self._calculate_ip_and_ea else ""

        cmd = (
            f"{memory} {self._xtb_path} "
            f"{xyz} --gfn {self._gfn_version} "
            f"{hess} {vipea} --parallel {self._num_cores} "
            f"--etemp {self._electronic_temperature} "
            f"{solvent} --chrg {self._charge} "
            f"--uhf {self._num_unpaired_electrons} -I det_control.in"
        )

        try:
            os.chdir(output_dir)
            self._write_detailed_control()
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
        xyz = output_dir / "input_structure.xyz"
        out_file = output_dir / "energy.output"
        mol.write(xyz)
        self._run_xtb(
            xyz=xyz,
            out_file=out_file,
            init_dir=init_dir,
            output_dir=output_dir,
        )
        yield

    def get_results(self, mol: stk.Molecule) -> XTBResults:
        """Calculate the xTB properties of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The properties, with units, from xTB calculations.

        """
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        out_file = output_dir / "energy.output"

        return XTBResults(
            generator=self.calculate(mol),
            output_file=out_file,
        )

    def get_energy(self, mol: stk.Molecule) -> float:
        """Calculate the energy of `mol`.

        Parameters:
            mol:
                The :class:`stk.Molecule` whose energy is to be calculated.

        Returns:
            The energy.

        """
        return self.get_results(mol).get_total_energy()[0]
