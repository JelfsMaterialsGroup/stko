import logging
import os
import shutil
import subprocess as sp
import uuid
from collections.abc import Iterable
from pathlib import Path

import stk

from stko._internal.calculators.extractors.xtb_extractor import XTBExtractor
from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT
from stko._internal.utilities.exceptions import (
    ConvergenceError,
    InvalidSolventError,
    NotCompletedError,
    NotStartedError,
    OptimizerError,
    PathError,
    SettingConflictError,
)
from stko._internal.utilities.utilities import is_valid_xtb_solvent

logger = logging.getLogger(__name__)


class XTB(Optimizer):
    """Uses GFN-xTB to optimize molecules.

    See Also:
        * https://xtb-docs.readthedocs.io/en/latest/setup.html

    Notes:
        When running :meth:`optimize`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        Furthermore, :meth:`optimize` will check that the
        structure is adequately optimized by checking for negative
        frequencies after a Hessian calculation. `max_runs` can be
        provided to the initializer to set the maximum number of
        optimizations which will be attempted at the given
        `opt_level` to obtain an optimized structure. However, we outline
        in the examples how to iterate over `opt_levels` to increase
        convergence criteria and hopefully obtain an optimized structure.
        The presence of negative frequencies can occur even when the
        optimization has converged based on the given `opt_level`.

        We thank Andrew Tarzia and Alejandro Santana-Bonilla for their
        contributions to this code.

    Parameters:
        xtb_path:
            The path to the xTB executable.

        gfn_version:
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir:
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level:
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        max_runs:
            Maximum number of optimizations to attempt in a row.

        calculate_hessian:
            Toggle calculation of the hessian and vibrational
            frequencies after optimization. ``True`` is required to
            check that the structure is completely optimized.
            ``False`` will drastically speed up the calculation but
            potentially provide incomplete optimizations and forces
            :attr:`max_runs` to be ``1``.

        num_cores:
            The number of cores xTB should use.

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
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        write_sasa_info:
            If ``True``, the detailed info input will request ``gbsa=True`` and
            output SASA information from xtb. Requires a solvent model to be
            used.

    Examples:
        Note that for :class:`.ConstructedMolecule` objects constructed by
        ``stk``, :class:`XTB` should usually be used in a
        :class:`.OptimizerSequence`. This is because xTB only uses
        xyz coordinates as input and so will not recognize the long bonds
        created during construction. An optimizer which can minimize
        these bonds should be used before :class:`XTB`.

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

            xtb = stko.OptimizerSequence(
                stko.UFF(),
                stko.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
            )
            polymer = xtb.optimize(polymer)


        By default, all optimizations with xTB are performed using the
        ``--ohess`` flag, which forces the calculation of a numerical
        Hessian, thermodynamic properties and vibrational frequencies.
        :meth:`optimize` will check that the structure is appropriately
        optimized (i.e. convergence is obtained and no negative vibrational
        frequencies are present) and continue optimizing a structure (up to
        `max_runs` times) until this is achieved. This loop, by
        default, will be performed at the same `opt_level`. The
        following example shows how a user may optimize structures with
        tigher convergence criteria (i.e. different `opt_level`)
        until the structure is sufficiently optimized. Furthermore, the
        calculation of the Hessian can be turned off using
        `max_runs` to ``1`` and `calculate_hessian` to ``False``.

        .. code-block:: python

            # Use crude optimization with max_runs=1 because this will
            # not achieve optimization and rerunning it is unproductive.
            xtb_crude = stko.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                output_dir='xtb_crude',
                unlimited_memory=True,
                opt_level='crude',
                max_runs=1,
                calculate_hessian=True
            )
            # Use normal optimization with max_runs == 2.
            xtb_normal = stko.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                output_dir='xtb_normal',
                unlimited_memory=True,
                opt_level='normal',
                max_runs=2
            )
            # Use vtight optimization with max_runs == 2, which should
            # achieve sufficient optimization.
            xtb_vtight = stko.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                output_dir='xtb_vtight',
                unlimited_memory=True,
                opt_level='vtight',
                max_runs=2
            )

            optimizers = [xtb_crude, xtb_normal, xtb_vtight]
            for optimizer in optimizers:
                polymer = optimizer.optimize(polymer)
                if polymer not in optimizer.incomplete:
                    break

    """

    incomplete: set[stk.Molecule]
    """
    Molecules passed to :meth:`optimize` whose optimzation was incomplete.
    """

    def __init__(  # noqa: PLR0913
        self,
        xtb_path: Path | str,
        gfn_version: int = 2,
        output_dir: Path | str | None = None,
        opt_level: str = "normal",
        max_runs: int = 2,
        calculate_hessian: bool = True,
        num_cores: int = 1,
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
                msg = (
                    "XTB: No solvent valid for version",
                    f" {gfn_version!r}.",
                )
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

        if not calculate_hessian and max_runs != 1:
            max_runs = 1
            logger.warning(
                "Requested that hessian calculation was skipped "
                "but the number of optimizations requested was "
                "greater than 1. The number of optimizations has been "
                "set to 1."
            )

        self._check_path(xtb_path)
        self._xtb_path = Path(xtb_path)
        self._gfn_version = str(gfn_version)
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._opt_level = opt_level
        self._max_runs = max_runs
        self._calculate_hessian = calculate_hessian
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_model = solvent_model
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self._write_sasa_info = write_sasa_info
        self.incomplete = set()

    def _check_path(self, path: Path | str) -> None:
        path = Path(path)
        if not path.exists():
            msg = f"XTB not found at {path}"
            raise PathError(msg)

    def _has_neg_frequencies(self, output_file: Path | str) -> bool:
        """Check for negative frequencies.

        Parameters:
            output_file:
                Name of output file with xTB results.

        Returns:
            ``True`` if a negative frequency is present.

        """
        xtbext = XTBExtractor(output_file=output_file)
        # Check for one negative frequency, excluding the first
        # 6 frequencies.
        return any(x < 0 for x in xtbext.frequencies[6:])

    def _is_complete(self, output_file: Path | str) -> bool:
        """Check if xTB optimization has completed and converged.

        Parameters:
            output_file:
                Name of xTB output file.

        Returns:
            ``False`` if a negative frequency is present.

        Raises:
            :class:`XTBOptimizerError`:
                if the optimization failed.

            :class:`XTBConvergenceError`:
                if the optimization did not converge.

        """
        output_file = Path(output_file)
        if not output_file.exists():
            # No simulation has been run.
            msg = "XTB: Optimization failed to start"
            raise OptimizerError(msg)

        # If convergence is achieved, then .xtboptok should exist.
        if Path(".xtboptok").exists():
            # Check for negative frequencies in output file if the
            # hessian was calculated.
            # Return True if there exists at least one.
            if self._calculate_hessian:
                return not self._has_neg_frequencies(output_file)
            return True

        if Path("NOT_CONVERGED").exists():
            msg = "XTB: Optimization not converged."
            raise ConvergenceError(msg)
        msg = "XTB: Optimization failed to complete"
        raise OptimizerError(msg)

    def _run_xtb(self, xyz: str, out_file: Path | str) -> None:
        """Run GFN-xTB.

        Parameters
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

        """
        out_file = Path(out_file)

        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        # Set optimization level and type.
        if self._calculate_hessian:
            # Do optimization and check hessian.
            optimization = f"--ohess {self._opt_level}"
        else:
            # Do optimization.
            optimization = f"--opt {self._opt_level}"

        if self._solvent is not None:
            solvent = f"--{self._solvent_model} {self._solvent} "
        else:
            solvent = ""

        cmd = (
            f"{memory} {self._xtb_path} {xyz} "
            f"--gfn {self._gfn_version} "
            f"{optimization} --parallel {self._num_cores} "
            f"--etemp {self._electronic_temperature} "
            f"{solvent} --chrg {self._charge} "
            f"--uhf {self._num_unpaired_electrons} -I det_control.in"
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _write_detailed_control(self) -> None:
        sasa_info = "$write\n gbsa=true\n" if self._write_sasa_info else ""
        string = f"$gbsa\n gbsagrid={self._solvent_grid}\n{sasa_info}"

        with Path("det_control.in").open("w") as f:
            f.write(string)

    def _run_optimizations(
        self,
        mol: MoleculeT,
    ) -> tuple[MoleculeT, bool]:
        """Run loop of optimizations on `mol` using xTB.

        Parameters
            mol:
                The molecule to be optimized.

        Returns:
            The optimized molecule and ``True`` if the calculation
            is complete or ``False`` if the calculation is incomplete.

        """
        for run in range(self._max_runs):
            xyz = f"input_structure_{run+1}.xyz"
            out_file = f"optimization_{run+1}.output"
            mol.write(xyz)
            self._write_detailed_control()
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = Path("xtbhess.coord")
            output_xyz = Path("xtbopt.xyz")
            opt_complete = self._is_complete(out_file)
            if not opt_complete:
                if coord_file.exists():
                    # The calculation is incomplete.
                    # Update mol from xtbhess.coord and continue.
                    mol = mol.with_structure_from_file(coord_file)
                else:
                    # Update mol from xtbopt.xyz.
                    mol = mol.with_structure_from_file(output_xyz)
                    # If the negative frequencies are small, then GFN
                    # may not produce the restart file. If that is the
                    # case, exit optimization loop and warn.
                    self.incomplete.add(mol)
                    msg = f"Small negative frequencies present in {mol}."
                    logging.warning(msg)
                    return mol, opt_complete
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol = mol.with_structure_from_file(output_xyz)
                break

        return mol, opt_complete

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol:
                The molecule to be optimized.

        Returns:
            The optimized molecule.

        """
        # Remove mol from self.incomplete if present.
        if mol in self.incomplete:
            self.incomplete.remove(mol)

        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        if output_dir.exists():
            shutil.rmtree(output_dir)

        output_dir.mkdir(parents=True)
        init_dir = Path.cwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimizations(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(mol)
            msg = f"Optimization is incomplete for {mol}."
            logging.warning(msg)

        return mol


class XTBCREST(Optimizer):
    """Uses GFN-n to run CREST on molecules.

    See Also:
        * GFN-n: https://xtb-docs.readthedocs.io/en/latest/setup.html
        * CREST: https://xtb-docs.readthedocs.io/en/latest/crestcmd.html

    Parameters:
        crest_path:
            The path to the CREST executable.

        xtb_path:
            The path to the xTB executable.
            Version >6.3.0 is required.

        gfn_version:
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir:
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level:
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        md_len:
            Set length of the meta-dynamics simulations (MTD) in ps.
            Default is chosen based on size and flexibility of the
            system.

        ewin:
            Set the energy threshold in kcal/mol for conformer
            selection. Double this is used in crude optimization.
            Defaults ot 5 kcal/mol and is overridden by
            :attr:`speed_setting`.

        speed_setting:
            Conformer search speed setting. Fast methods turn off
            parts of the calculations and alter MD run times.
            Defaults to no modification of iMTD-GC algorithm: `None`.
            Can be one of ``'norotmd'``, ``'quick'``, ``'squick'``
            or ``'mquick'``.
            Overrides :attr:`ewin` with chosen parameters.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        keepdir:
            `True` to keep subdirectories from MD runs.
            Defaults to `False`.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        num_cores:
            The number of cores CREST should use.

        charge:
            Formal molecular charge.

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

        num_unpaired_electrons:
            Number of unpaired electrons.

        cross:
            Whether or not structure crossing is performed.

        unlimited_memory:
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

    Notes:
        Requires version > 6.2 of xtb.

        When running :meth:`optimize`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        *Restrictions to iMTD-GC Algorithm*
        Z-matrix sorting is forced to be off because :class:`stk.Molecules`
        cannot have their atom ordering changed by an external program at
        this stage.

    Examples:
        Note that for :class:`.ConstructedMolecule` objects constructed by
        ``stk``, :class:`XTBCREST` should usually be used in a
        :class:`.OptimizerSequence`. This is because xTB only uses
        xyz coordinates as input and so will not recognize the long bonds
        created during construction. An optimizer which can minimize
        these bonds should be used before :class:`XTBCREST`. Further,
        CREST runs best on an input structure optimized at the same level
        used throughout the algorithm (i.e. :class:`XTB`).

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

            xtb = stko.OptimizerSequence(
                stko.UFF(),
                stko.XTB(
                    xtb_path='/opt/gfnxtb/xtb',
                    unlimited_memory=True,
                ),
                # Perform quick conformer search.
                stko.XTBCREST(
                    crest_path='/opt/crest/crest',
                    xtb_path='/opt/gfnxtb/xtb',
                    opt_level='normal',
                    speed_setting='quick',
                    unlimited_memory=True,
                ),
            )
            polymer = xtb.optimize(polymer)

        While most auxiliary and MD run files will be deleted if
        :attr:`keepdir` is `False`, the files listed below are kept in
        :attr:`output_dir` and could be useful for further analysis!

            | crest_best.xyz: Structure of the lowest energy conformer.
            |    This structure is output by :meth:`optimize`.
            | crest_conformers.xyz: Structure of all conformers.
            |    All conformers with RMSD and energy threshold of the
            |    lowest energy conformer.
            | crest_rotamers.xyz: Structure of all rotamers.
            |    All unique rotamers explored by CREST.
            | crest.energies: Relative conformer energies in a.u..
            | gfn_topo: GFN-FF binary topology file.
            |    Defines the molecules force field topology.


    """

    def __init__(  # noqa: PLR0913
        self,
        crest_path: str,
        xtb_path: str,
        gfn_version: int = 2,
        output_dir: str | None = None,
        opt_level: str = "normal",
        md_len: float | None = None,
        ewin: float = 5,
        speed_setting: None | str = None,
        keepdir: bool = False,
        num_cores: int = 4,
        cross: bool = True,
        charge: int = 0,
        electronic_temperature: float = 300,
        solvent_model: str = "gbsa",
        solvent: str | None = None,
        num_unpaired_electrons: int = 0,
        unlimited_memory: bool = False,
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

        self._check_path(crest_path)
        self._check_path(xtb_path)
        self._crest_path = crest_path
        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._opt_level = opt_level
        self._mdlen = md_len

        if speed_setting is not None and ewin != 5:  # noqa: PLR2004
            msg = (
                f"CREST: The chosen speed setting {speed_setting} will ",
                f"override the chosen energy window {ewin}.",
            )
            raise SettingConflictError(msg)

        self._ewin = ewin
        self._speed_setting = speed_setting
        self._keepdir = keepdir
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent_model = solvent_model
        self._solvent = solvent
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._cross = cross
        self._unlimited_memory = unlimited_memory

    def _check_path(self, path: Path | str) -> None:
        path = Path(path)
        if not path.exists():
            msg = f"XTB or CREST not found at {path}"
            raise PathError(msg)

    def _is_complete(
        self, output_file: Path | str, output_xyzs: Iterable[Path]
    ) -> bool:
        """Check if CREST run has completed.

        Parameters:
            output_file:
                Name of CREST output file.

            output_xyzs:
                Name of CREST conformer output files.
                crest_best.xyz > Best conformer, exists throughout run.
                crest_conformers.xyz > All conformers,
                    exists throughout run.

        Returns:
            ``False`` if a negative frequency is present.

        Raises:
            :class:`CRESTNotStartedError`:
                if the CREST run failed to start.

            :class:`CRESTNotCompletedError`:
                if the CREST run failed to complete.

        """
        output_file = Path(output_file)
        if not output_file.exists():
            # No simulation has been run.
            msg = "CREST run did not start"
            raise NotStartedError(msg)

        if any(not xyz.exists() for xyz in output_xyzs):
            # Best conformer was not output.
            msg = "CREST run did not complete"
            raise NotCompletedError(msg)

        return True

    def _run_crest(self, xyz: str, out_file: Path | str) -> None:
        """Run CREST along side GFN-xTB.

        Parameters:
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

        """
        out_file = Path(out_file)

        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        if self._solvent is not None:
            solvent = f"--{self._solvent_model} {self._solvent}"
        else:
            solvent = ""

        # Set optimization level and type.
        optimization = f"-opt {self._opt_level}"
        mdlen = "" if self._mdlen is None else f"-mdlen {self._mdlen} "
        keepdirs = "-keepdir" if self._keepdir is True else ""
        if self._speed_setting is not None:
            speed_settings = f"-{self._speed_setting}"
            ewin = ""
        else:
            speed_settings = ""
            ewin = f"-ewin {self._ewin} "
        cross = "-nocross" if self._cross is False else ""

        cmd = (
            f"{memory} {self._crest_path} {xyz} "
            f"-xnam {self._xtb_path} -nozs {cross} "
            f"{ewin} {mdlen}"
            f"{solvent} -chrg {self._charge} "
            f"-uhf {self._num_unpaired_electrons} "
            f"-gfn{self._gfn_version} "
            f"{speed_settings} "
            f"{optimization} -T {self._num_cores} {keepdirs}"
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _run_optimization(
        self,
        mol: MoleculeT,
    ) -> tuple[MoleculeT, bool]:
        """Run loop of optimizations on `mol` using xTB.

        Parameters:
            mol:
                The molecule to be optimized.

        Returns:
            The optimized molecule and ``True`` if the calculation
            is complete or ``False`` if the calculation is incomplete.

        """
        xyz = "input_structure.xyz"
        out_file = "crest.output"
        mol.write(xyz)
        self._run_crest(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyzs = [
            Path("crest_best.xyz"),
            Path("crest_conformers.xyz"),
        ]
        opt_complete = self._is_complete(out_file, output_xyzs)
        mol = mol.with_structure_from_file(output_xyzs[0])

        return mol, opt_complete

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol:
                The molecule to be optimized.

        Returns:
            The optimized molecule.

        """
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        if output_dir.exists():
            shutil.rmtree(output_dir)

        output_dir.mkdir(parents=True)
        init_dir = Path.cwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            msg = f"CREST run is incomplete for {mol}."
            logging.warning(msg)

        return mol


class XTBFF(Optimizer):
    """Uses GFN-FF to optimize molecules.

    See Also:
        * GFN-FF: https://xtb-docs.readthedocs.io/en/latest/gfnff.html

    Parameters:
        xtb_path:
            The path to the xTB executable.

        output_dir:
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level:
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        num_cores:
            The number of cores xTB should use.

        charge:
            Formal molecular charge.

        unlimited_memory:
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.


    Notes:
        GFN-FF requires >= version 6.3 of xtb.

        When running :meth:`optimize`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        Currently, we only provide inputs that work with GFN-FF,
        specifically the charge of the system. Other electronic properties
        of the molecule are not relavent to a forcefield optimisation.

    Examples:
        Note that for :class:`.ConstructedMolecule` objects constructed by
        ``stk``, :class:`XTBFF` should usually be used in a
        :class:`.OptimizerSequence`. This is because xTB only uses
        xyz coordinates as input and so will not recognize the long bonds
        created during construction. An optimizer which can minimize
        these bonds should be used before :class:`XTBFF`.

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

            xtb = stko.OptimizerSequence(
                stko.UFF(),
                stko.XTBFF(
                    xtb_path='/opt/gfnxtb/xtb',
                    unlimited_memory=True,
                ),
            )
            polymer = xtb.optimize(polymer)


    """

    def __init__(  # noqa: PLR0913
        self,
        xtb_path: str,
        output_dir: Path | str | None = None,
        opt_level: str = "normal",
        num_cores: int = 1,
        charge: int = 0,
        unlimited_memory: bool = False,
    ) -> None:
        self._check_path(xtb_path)
        self._xtb_path = xtb_path
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._opt_level = opt_level
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._unlimited_memory = unlimited_memory

    def _check_path(self, path: Path | str) -> None:
        path = Path(path)
        if not path.exists():
            msg = f"XTB not found at {path}"
            raise PathError(msg)

    def _is_complete(self, output_file: Path | str) -> bool:
        """Check if xTB optimization has completed and converged.

        Parameters:
            output_file:
                Name of xTB output file.

        Returns:
            Returns ``False`` if a negative frequency is present.

        Raises:
            :class:`XTBOptimizerError`:
                if the optimization failed.

            :class:`XTBConvergenceError`:
                if the optimization did not converge.

        """
        output_file = Path(output_file)

        if not output_file.exists():
            # No simulation has been run.
            msg = "XTB: Optimization failed to start"
            raise OptimizerError(msg)
        # If convergence is achieved, then .xtboptok should exist.
        if Path(".xtboptok").exists():
            return True
        if Path("NOT_CONVERGED").exists():
            msg = "XTB: Optimization not converged."
            raise ConvergenceError(msg)
        msg = "XTB: Optimization failed to complete"
        raise OptimizerError(msg)

    def _run_xtb(self, xyz: str, out_file: Path | str) -> None:
        """Run GFN-xTB.

        Parameters:
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.


        """
        out_file = Path(out_file)

        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        # Set optimization level and type.
        optimization = f"--opt {self._opt_level}"

        cmd = (
            f"{memory} {self._xtb_path} {xyz} "
            f"--gfnff "
            f"{optimization} --parallel {self._num_cores} "
            f"--chrg {self._charge} "
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _run_optimization(
        self,
        mol: MoleculeT,
    ) -> tuple[MoleculeT, bool]:
        """Run loop of optimizations on `mol` using xTB.

        Parameters:
            mol: The molecule to be optimized.

        Returns:
            The optimized molecule and ``True`` if the calculation
            is complete or ``False`` if the calculation is incomplete.

        """
        xyz = "input_structure_ff.xyz"
        out_file = "optimization_ff.output"
        mol.write(xyz)
        self._run_xtb(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyz = "xtbopt.xyz"
        opt_complete = self._is_complete(out_file)
        mol = mol.with_structure_from_file(output_xyz)

        return mol, opt_complete

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol: The molecule to be optimized.

        Returns:
            The optimized molecule.

        """
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        if output_dir.exists():
            shutil.rmtree(output_dir)

        output_dir.mkdir(parents=True)
        init_dir = Path.cwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            msg = f"Optimization is incomplete for {mol}."
            logging.warning(msg)

        return mol


class XTBFFCREST(Optimizer):
    """Uses GFN-FF to run CREST on molecules.

    See Also:
        * GFN-FF: https://xtb-docs.readthedocs.io/en/latest/gfnff.html
        * CREST: https://xtb-docs.readthedocs.io/en/latest/crestcmd.html

    Parameters:

        crest_path:
            The path to the CREST executable.

        xtb_path:
            The path to the xTB executable.
            Version >6.3.0 is required.

        output_dir:
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level:
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        md_len:
            Set length of the meta-dynamics simulations (MTD) in ps.
            Default is chosen based on size and flexibility of the
            system.

        ewin:
            Set the energy threshold in kcal/mol for conformer
            selection. Double this is used in crude optimization.
            Defaults ot 5 kcal/mol and is overridden by
            :attr:`speed_setting`.

        speed_setting:
            Conformer search speed setting. Fast methods turn off
            parts of the calculations and alter MD run times.
            Defaults to no modification of iMTD-GC algorithm: `None`.
            Can be one of ``'norotmd'``, ``'quick'``, ``'squick'``
            or ``'mquick'``.
            Overrides :attr:`ewin` with chosen parameters.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        keepdir:
            `True` to keep subdirectories from MD runs.
            Defaults to `False`.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        num_cores:
            The number of cores CREST should use.

        charge:
            Formal molecular charge.

        cross:
            Whether or not structure crossing is performed.

        unlimited_memory:
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.


    Notes:
        GFN-FF requires version 6.3 of xtb.

        When running :meth:`optimize`, this calculator changes the
        present working directory with :func:`os.chdir`. The original
        working directory will be restored even if an error is raised, so
        unless multi-threading is being used this implementation detail
        should not matter.

        Currently, we only provide inputs that work with GFN-FF,
        specifically the charge of the system. Other electronic properties
        of the molecule are not relavent to a forcefield optimisation.
        We intend on adding more options in the future!

        *Restrictions to iMTD-GC Algorithm*
        Z-matrix sorting is forced to be off because :class:`stk.Molecules`
        cannot have their atom ordering changed by an external program at
        this stage.

    Examples:
        Note that for :class:`.ConstructedMolecule` objects constructed by
        ``stk``, :class:`XTBFFCREST` should usually be used in a
        :class:`.OptimizerSequence`. This is because xTB only uses
        xyz coordinates as input and so will not recognize the long bonds
        created during construction. An optimizer which can minimize
        these bonds should be used before :class:`XTBFFCREST`. Further,
        CREST runs best on an input structure optimized at the same level
        used throughout the algorithm (i.e. :class:`XTBFF`).

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

            xtb = stko.OptimizerSequence(
                stko.UFF(),
                stko.XTBFF(
                    xtb_path='/opt/gfnxtb/xtb',
                    unlimited_memory=True,
                ),
                # Perform quick conformer search.
                stko.XTBFFCREST(
                    crest_path='/opt/crest/crest',
                    xtb_path='/opt/gfnxtb/xtb',
                    opt_level='normal',
                    speed_setting='quick',
                    unlimited_memory=True,
                ),
            )
            polymer = xtb.optimize(polymer)

        While most auxiliary and MD run files will be deleted if
        :attr:`keepdir` is `False`, the files listed below are kept in
        :attr:`output_dir` and could be useful for further analysis!

            | crest_best.xyz: Structure of the lowest energy conformer.
            |    This structure is output by :meth:`optimize`.
            | crest_conformers.xyz: Structure of all conformers.
            |    All conformers with RMSD and energy threshold of the
            |    lowest energy conformer.
            | crest_rotamers.xyz: Structure of all rotamers.
            |    All unique rotamers explored by CREST.
            | crest.energies: Relative conformer energies in a.u..
            | gfn_topo: GFN-FF binary topology file.
            |    Defines the molecules force field topology.


    """

    def __init__(  # noqa: PLR0913
        self,
        crest_path: str,
        xtb_path: str,
        output_dir: Path | str | None = None,
        opt_level: str = "normal",
        md_len: float | None = None,
        ewin: float = 5,
        speed_setting: str | None = None,
        keepdir: bool = False,
        num_cores: int = 1,
        charge: int = 0,
        cross: bool = True,
        unlimited_memory: bool = False,
    ) -> None:
        self._check_path(crest_path)
        self._check_path(xtb_path)
        self._crest_path = crest_path
        self._xtb_path = xtb_path
        self._output_dir = None if output_dir is None else Path(output_dir)
        self._opt_level = opt_level
        self._mdlen = md_len

        if speed_setting is not None and ewin != 5:  # noqa: PLR2004
            msg = (
                f"CREST: The chosen speed setting {speed_setting} will ",
                f"override the chosen energy window {ewin}.",
            )
            raise SettingConflictError(msg)

        self._ewin = ewin
        self._speed_setting = speed_setting
        self._keepdir = keepdir
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._cross = cross
        self._unlimited_memory = unlimited_memory

    def _check_path(self, path: Path | str) -> None:
        path = Path(path)
        if not path.exists():
            msg = f"XTB or CREST not found at {path}"
            raise PathError(msg)

    def _is_complete(
        self, output_file: Path, output_xyzs: Iterable[Path]
    ) -> bool:
        """Check if CREST run has completed.

        Parameters:
            output_file:
                Name of CREST output file.

            output_xyzs:
                Name of CREST conformer output files.
                crest_best.xyz > Best conformer, exists throughout run.
                crest_conformers.xyz > All conformers,
                    exists throughout run.

        Returns:
            Returns ``False`` if a negative frequency is present.

        Raises:
            :class:`CRESTNotStartedError`:
                if the CREST run failed to start.

            :class:`CRESTNotCompletedError`:
                if the CREST run failed to complete.

        """
        if not output_file.exists():
            # No simulation has been run.
            msg = "CREST run did not start"
            raise NotStartedError(msg)

        if any(not xyz.exists() for xyz in output_xyzs):
            # Best conformer was not output.
            msg = "CREST run did not complete"
            raise NotCompletedError(msg)

        return True

    def _run_crest(self, xyz: Path, out_file: Path) -> None:
        """Run CREST along side GFN-xTB.

        Parameters:
            xyz:
                The name of the input structure ``.xyz`` file.

            out_file:
                The name of output file with xTB results.

        """
        # Modify the memory limit.
        memory = "ulimit -s unlimited ;" if self._unlimited_memory else ""

        # Set optimization level and type.
        optimization = f"-opt {self._opt_level}"
        mdlen = "" if self._mdlen is None else f"-mdlen {self._mdlen} "
        keepdirs = "-keepdir" if self._keepdir is True else ""
        if self._speed_setting is not None:
            speed_settings = f"-{self._speed_setting}"
            ewin = ""
        else:
            speed_settings = ""
            ewin = f"-ewin {self._ewin} "
        cross = "-nocross" if self._cross is False else ""

        cmd = (
            f"{memory} {self._crest_path} {xyz} "
            f"-xnam {self._xtb_path} -nozs {cross} "
            f"{ewin} {mdlen}"
            f"-chrg {self._charge} "
            f"-gff {speed_settings} "
            f"{optimization} -T {self._num_cores} {keepdirs}"
        )

        with out_file.open("w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(  # noqa: S602
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def _run_optimization(
        self,
        mol: MoleculeT,
    ) -> tuple[MoleculeT, bool]:
        """Run loop of optimizations on `mol` using xTB.

        Parameters:
            mol:
                The molecule to be optimized.

        Returns:
            mol:
                The optimized molecule.

            opt_complete:
                Returns ``True`` if the calculation is complete and
                ``False`` if the calculation is incomplete.

        """
        xyz = Path("input_structure_ff.xyz")
        out_file = Path("crest_ff.output")
        mol.write(xyz)
        self._run_crest(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyzs = [
            Path("crest_best.xyz"),
            Path("crest_conformers.xyz"),
        ]
        opt_complete = self._is_complete(out_file, output_xyzs)
        mol = mol.with_structure_from_file(output_xyzs[0])

        return mol, opt_complete

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol:
                The molecule to be optimized.

        Returns:
            The optimized molecule.

        """
        if self._output_dir is None:
            output_dir = Path(str(uuid.uuid4().int)).resolve()
        else:
            output_dir = self._output_dir.resolve()

        if output_dir.exists():
            shutil.rmtree(output_dir)

        output_dir.mkdir(parents=True)
        init_dir = Path.cwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            msg = f"CREST run is incomplete for {mol}."
            logging.warning(msg)

        return mol
