"""
XTB Optimizers
==============

#. :class:`.XTB`
#. :class:`.XTBCREST`
#. :class:`.XTBFF`
#. :class:`.XTBFFCREST`

Wrappers for optimizers within the :mod:`xtb` code.

"""

import logging
import os
import shutil
import uuid
import subprocess as sp

from .optimizers import Optimizer
from ..utilities import (
    is_valid_xtb_solvent,
    XTBInvalidSolventError,
    XTBExtractor
)

logger = logging.getLogger(__name__)


class XTBOptimizerError(Exception):
    ...


class XTBConvergenceError(XTBOptimizerError):
    ...


class CRESTOptimizerError(Exception):
    ...


class CRESTNotStartedError(CRESTOptimizerError):
    ...


class CRESTNotCompletedError(CRESTOptimizerError):
    ...


class CRESTSettingConflictError(CRESTOptimizerError):
    ...


class XTB(Optimizer):
    """
    Uses GFN-xTB [1]_ to optimize molecules.

    Notes
    -----
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

    *Contributors*
    We thank Andrew Tarzia and Alejandro Santana-Bonilla for their
    contributions to this code.

    Attributes
    ----------
    incomplete : :class:`set` of :class:`.Molecule`
        A :class:`set` of molecules passed to :meth:`optimize` whose
        optimzation was incomplete.

    Examples
    --------
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

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(
        self,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        opt_level='normal',
        max_runs=2,
        calculate_hessian=True,
        num_cores=1,
        electronic_temperature=300,
        solvent_model='gbsa',
        solvent=None,
        solvent_grid='normal',
        charge=0,
        num_unpaired_electrons=0,
        unlimited_memory=False,
    ):
        """
        Initialize a :class:`XTB` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`int`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        max_runs : :class:`int`, optional
            Maximum number of optimizations to attempt in a row.

        calculate_hessian : :class:`bool`, optional
            Toggle calculation of the hessian and vibrational
            frequencies after optimization. ``True`` is required to
            check that the structure is completely optimized.
            ``False`` will drastically speed up the calculation but
            potentially provide incomplete optimizations and forces
            :attr:`max_runs` to be ``1``.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

        solvent_model : :class:`str`
            Solvent model to use out of older `gbsa` and newer `alpb`.
            `gbsa` is default for backwards compatability, but `alpb`
            is recommended.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        charge : :class:`int`, optional
            Formal molecular charge.

        num_unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        """

        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(
                gfn_version=gfn_version,
                solvent_model=solvent_model,
                solvent=solvent,
            ):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} and model {solvent_model!r}',
                    f' is invalid for version {gfn_version!r}.'
                )

        if not calculate_hessian and max_runs != 1:
            max_runs = 1
            logger.warning(
                'Requested that hessian calculation was skipped '
                'but the number of optimizations requested was '
                'greater than 1. The number of optimizations has been '
                'set to 1.'
            )

        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
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
        self.incomplete = set()

    def _has_neg_frequencies(self, output_file):
        """
        Check for negative frequencies.

        Parameters
        ----------
        output_file : :class:`str`
            Name of output file with xTB results.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if a negative frequency is present.

        """
        xtbext = XTBExtractor(output_file=output_file)
        # Check for one negative frequency, excluding the first
        # 6 frequencies.
        return any(x < 0 for x in xtbext.frequencies[6:])

    def _is_complete(self, output_file):
        """
        Check if xTB optimization has completed and converged.

        Parameters
        ----------
        output_file : :class:`str`
            Name of xTB output file.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`XTBOptimizerError`
            If the optimization failed.

        :class:`XTBConvergenceError`
            If the optimization did not converge.

        """
        if not os.path.exists(output_file):
            # No simulation has been run.
            raise XTBOptimizerError('Optimization failed to start')

        # If convergence is achieved, then .xtboptok should exist.
        if os.path.exists('.xtboptok'):
            # Check for negative frequencies in output file if the
            # hessian was calculated.
            # Return True if there exists at least one.
            if self._calculate_hessian:
                return not self._has_neg_frequencies(output_file)
            else:
                return True

        elif os.path.exists('NOT_CONVERGED'):
            raise XTBConvergenceError('Optimization not converged.')
        else:
            raise XTBOptimizerError('Optimization failed to complete')

    def _run_xtb(self, xyz, out_file):
        """
        Run GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        if self._calculate_hessian:
            # Do optimization and check hessian.
            optimization = f'--ohess {self._opt_level}'
        else:
            # Do optimization.
            optimization = f'--opt {self._opt_level}'

        if self._solvent is not None:
            solvent = f'--{self._solvent_model} {self._solvent} '
        else:
            solvent = ''

        cmd = (
            f'{memory} {self._xtb_path} {xyz} '
            f'--gfn {self._gfn_version} '
            f'{optimization} --parallel {self._num_cores} '
            f'--etemp {self._electronic_temperature} '
            f'{solvent} --chrg {self._charge} '
            f'--uhf {self._num_unpaired_electrons} -I det_control.in'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _write_detailed_control(self):
        string = f'$gbsa\n   gbsagrid={self._solvent_grid}'

        with open('det_control.in', 'w') as f:
            f.write(string)

    def _run_optimizations(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        opt_complete : :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """

        for run in range(self._max_runs):
            xyz = f'input_structure_{run+1}.xyz'
            out_file = f'optimization_{run+1}.output'
            mol.write(xyz)
            self._write_detailed_control()
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = 'xtbhess.coord'
            output_xyz = 'xtbopt.xyz'
            opt_complete = self._is_complete(out_file)
            if not opt_complete:
                if os.path.exists(coord_file):
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
                    logging.warning(
                        f'Small negative frequencies present in {mol}.'
                    )
                    return mol, opt_complete
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol = mol.with_structure_from_file(output_xyz)
                break

        return mol, opt_complete

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

        # Remove mol from self.incomplete if present.
        if mol in self.incomplete:
            self.incomplete.remove(mol)

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimizations(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(mol)
            logging.warning(f'Optimization is incomplete for {mol}.')

        return mol


class XTBCREST(Optimizer):
    """
    Uses GFN-n [1]_ to run CREST [2]_ on molecules.

    Notes
    -----
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

    Examples
    --------

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

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html
    .. [2] https://xtb-docs.readthedocs.io/en/latest/crestcmd.html

    """

    def __init__(
        self,
        crest_path,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        opt_level='normal',
        md_len=None,
        ewin=5,
        speed_setting=None,
        keepdir=False,
        num_cores=4,
        cross=True,
        charge=0,
        electronic_temperature=300,
        solvent_model='gbsa',
        solvent=None,
        num_unpaired_electrons=0,
        unlimited_memory=False,
    ):
        """
        Initialize a :class:`XTBCREST` instance.

        Parameters
        ----------
        crest_path : :class:`str`
            The path to the CREST executable.

        xtb_path : :class:`str`
            The path to the xTB executable.
            Version >6.3.0 is required.

        gfn_version : :class:`int`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        md_len : :class:`float`, optional.
            Set length of the meta-dynamics simulations (MTD) in ps.
            Default is chosen based on size and flexibility of the
            system.

        ewin : :class:`float`, optional.
            Set the energy threshold in kcal/mol for conformer
            selection. Double this is used in crude optimization.
            Defaults ot 5 kcal/mol and is overridden by
            :attr:`speed_setting`.

        speed_setting : :class:`str` or :class:`NoneType`, optional
            Conformer search speed setting. Fast methods turn off
            parts of the calculations and alter MD run times.
            Defaults to no modification of iMTD-GC algorithm: `None`.
            Can be one of ``'norotmd'``, ``'quick'``, ``'squick'``
            or ``'mquick'``.
            Overrides :attr:`ewin` with chosen parameters.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        keepdir : :class:`bool`, optional
            `True` to keep subdirectories from MD runs.
            Defaults to `False`.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        num_cores : :class:`int`, optional
            The number of cores CREST should use.

        charge : :class:`int`, optional
            Formal molecular charge.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

        solvent_model : :class:`str`
            Solvent model to use out of older `gbsa` and newer `alpb`.
            `gbsa` is default for backwards compatability, but `alpb`
            is recommended.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        num_unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        cross : :class:`bool`, optional
            Whether or not structure crossing is performed.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        """

        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(
                gfn_version=gfn_version,
                solvent_model=solvent_model,
                solvent=solvent,
            ):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} and model {solvent_model!r}',
                    f' is invalid for version {gfn_version!r}.'
                )

        self._crest_path = crest_path
        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
        self._opt_level = opt_level
        self._mdlen = md_len

        if speed_setting is not None and ewin != 5:
            raise CRESTSettingConflictError(
                f'The chosen speed setting {speed_setting} will ',
                f'override the chosen energy window {ewin}.'
            )

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

    def _is_complete(self, output_file, output_xyzs):
        """
        Check if CREST run has completed.

        Parameters
        ----------
        output_file : :class:`str`
            Name of CREST output file.

        output_xyzs : :class:`str`
            Name of CREST conformer output files.
            crest_best.xyz > Best conformer, exists throughout run.
            crest_conformers.xyz > All conformers,
                exists throughout run.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`CRESTNotStartedError`
            If the CREST run failed to start.

        :class:`CRESTNotCompletedError`
            If the CREST run failed to complete.

        """

        if not os.path.exists(output_file):
            # No simulation has been run.
            raise CRESTNotStartedError('CREST run did not start')

        elif any(not os.path.exists(i) for i in output_xyzs):
            # Best conformer was not output.
            raise CRESTNotCompletedError('CREST run did not complete')

        return True

    def _run_crest(self, xyz, out_file):
        """
        Run CREST along side GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        if self._solvent is not None:
            solvent = f'--{self._solvent_model} {self._solvent}'
        else:
            solvent = ''

        # Set optimization level and type.
        optimization = f'-opt {self._opt_level}'
        mdlen = '' if self._mdlen is None else f'-mdlen {self._mdlen} '
        keepdirs = '-keepdir' if self._keepdir is True else ''
        if self._speed_setting is not None:
            speed_settings = f'-{self._speed_setting}'
            ewin = ''
        else:
            speed_settings = ''
            ewin = f'-ewin {self._ewin} '
        cross = ('-nocross' if self._cross is False else '')

        cmd = (
            f'{memory} {self._crest_path} {xyz} '
            f'-xnam {self._xtb_path} -nozs {cross} '
            f'{ewin} {mdlen}'
            f'{solvent} -chrg {self._charge} '
            f'-uhf {self._num_unpaired_electrons} '
            f'-gfn{self._gfn_version} '
            f'{speed_settings} '
            f'{optimization} -T {self._num_cores} {keepdirs}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _run_optimization(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        opt_complete : :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """

        xyz = 'input_structure.xyz'
        out_file = 'crest.output'
        mol.write(xyz)
        self._run_crest(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyzs = [
            'crest_best.xyz',
            'crest_conformers.xyz',
        ]
        opt_complete = self._is_complete(out_file, output_xyzs)
        mol = mol.with_structure_from_file(output_xyzs[0])

        return mol, opt_complete

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

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            logging.warning(f'CREST run is incomplete for {mol}.')

        return mol


class XTBFF(Optimizer):
    """
    Uses GFN-FF [1]_ to optimize molecules.

    Notes
    -----
    GFN-FF requires >= version 6.3 of xtb.

    When running :meth:`optimize`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    Currently, we only provide inputs that work with GFN-FF,
    specifically the charge of the system. Other electronic properties
    of the molecule are not relavent to a forcefield optimisation.

    Examples
    --------
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

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/gfnff.html

    """

    def __init__(
        self,
        xtb_path,
        output_dir=None,
        opt_level='normal',
        num_cores=1,
        charge=0,
        unlimited_memory=False,
    ):
        """
        Initialize a :class:`XTB` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        charge : :class:`int`, optional
            Formal molecular charge.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        """

        self._xtb_path = xtb_path
        self._output_dir = output_dir
        self._opt_level = opt_level
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._unlimited_memory = unlimited_memory

    def _is_complete(self, output_file):
        """
        Check if xTB optimization has completed and converged.

        Parameters
        ----------
        output_file : :class:`str`
            Name of xTB output file.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`XTBOptimizerError`
            If the optimization failed.

        :class:`XTBConvergenceError`
            If the optimization did not converge.

        """
        if not os.path.exists(output_file):
            # No simulation has been run.
            raise XTBOptimizerError('Optimization failed to start')
        # If convergence is achieved, then .xtboptok should exist.
        if os.path.exists('.xtboptok'):
            return True
        elif os.path.exists('NOT_CONVERGED'):
            raise XTBConvergenceError('Optimization not converged.')
        else:
            raise XTBOptimizerError('Optimization failed to complete')

    def _run_xtb(self, xyz, out_file):
        """
        Run GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        optimization = f'--opt {self._opt_level}'

        cmd = (
            f'{memory} {self._xtb_path} {xyz} '
            f'--gfnff '
            f'{optimization} --parallel {self._num_cores} '
            f'--chrg {self._charge} '
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _run_optimization(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        opt_complete : :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """

        xyz = 'input_structure_ff.xyz'
        out_file = 'optimization_ff.output'
        mol.write(xyz)
        self._run_xtb(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyz = 'xtbopt.xyz'
        opt_complete = self._is_complete(out_file)
        mol = mol.with_structure_from_file(output_xyz)

        return mol, opt_complete

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

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            logging.warning(f'Optimization is incomplete for {mol}.')

        return mol


class XTBFFCREST(Optimizer):
    """
    Uses GFN-FF [1]_ to run CREST [2]_ on molecules.

    Notes
    -----
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

    Examples
    --------

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

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/gfnff.html
    .. [2] https://xtb-docs.readthedocs.io/en/latest/crestcmd.html

    """

    def __init__(
        self,
        crest_path,
        xtb_path,
        output_dir=None,
        opt_level='normal',
        md_len=None,
        ewin=5,
        speed_setting=None,
        keepdir=False,
        num_cores=1,
        charge=0,
        cross=True,
        unlimited_memory=False,
    ):
        """
        Initialize a :class:`XTBFFCREST` instance.

        Parameters
        ----------
        crest_path : :class:`str`
            The path to the CREST executable.

        xtb_path : :class:`str`
            The path to the xTB executable.
            Version >6.3.0 is required.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        md_len : :class:`float`, optional.
            Set length of the meta-dynamics simulations (MTD) in ps.
            Default is chosen based on size and flexibility of the
            system.

        ewin : :class:`float`, optional.
            Set the energy threshold in kcal/mol for conformer
            selection. Double this is used in crude optimization.
            Defaults ot 5 kcal/mol and is overridden by
            :attr:`speed_setting`.

        speed_setting : :class:`str` or :class:`NoneType`, optional
            Conformer search speed setting. Fast methods turn off
            parts of the calculations and alter MD run times.
            Defaults to no modification of iMTD-GC algorithm: `None`.
            Can be one of ``'norotmd'``, ``'quick'``, ``'squick'``
            or ``'mquick'``.
            Overrides :attr:`ewin` with chosen parameters.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        keepdir : :class:`bool`, optional
            `True` to keep subdirectories from MD runs.
            Defaults to `False`.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/crestcmd.html.

        num_cores : :class:`int`, optional
            The number of cores CREST should use.

        charge : :class:`int`, optional
            Formal molecular charge.

        cross : :class:`bool`, optional
            Whether or not structure crossing is performed.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        """

        self._crest_path = crest_path
        self._xtb_path = xtb_path
        self._output_dir = output_dir
        self._opt_level = opt_level
        self._mdlen = md_len

        if speed_setting is not None and ewin != 5:
            raise CRESTSettingConflictError(
                f'The chosen speed setting {speed_setting} will ',
                f'override the chosen energy window {ewin}.'
            )

        self._ewin = ewin
        self._speed_setting = speed_setting
        self._keepdir = keepdir
        self._num_cores = str(num_cores)
        self._charge = str(charge)
        self._cross = cross
        self._unlimited_memory = unlimited_memory

    def _is_complete(self, output_file, output_xyzs):
        """
        Check if CREST run has completed.

        Parameters
        ----------
        output_file : :class:`str`
            Name of CREST output file.

        output_xyzs : :class:`str`
            Name of CREST conformer output files.
            crest_best.xyz > Best conformer, exists throughout run.
            crest_conformers.xyz > All conformers,
                exists throughout run.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`CRESTNotStartedError`
            If the CREST run failed to start.

        :class:`CRESTNotCompletedError`
            If the CREST run failed to complete.

        """

        if not os.path.exists(output_file):
            # No simulation has been run.
            raise CRESTNotStartedError('CREST run did not start')

        elif any(not os.path.exists(i) for i in output_xyzs):
            # Best conformer was not output.
            raise CRESTNotCompletedError('CREST run did not complete')

        return True

    def _run_crest(self, xyz, out_file):
        """
        Run CREST along side GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        optimization = f'-opt {self._opt_level}'
        mdlen = '' if self._mdlen is None else f'-mdlen {self._mdlen} '
        keepdirs = '-keepdir' if self._keepdir is True else ''
        if self._speed_setting is not None:
            speed_settings = f'-{self._speed_setting}'
            ewin = ''
        else:
            speed_settings = ''
            ewin = f'-ewin {self._ewin} '
        cross = ('-nocross' if self._cross is False else '')

        cmd = (
            f'{memory} {self._crest_path} {xyz} '
            f'-xnam {self._xtb_path} -nozs {cross} '
            f'{ewin} {mdlen}'
            f'-chrg {self._charge} '
            f'-gff {speed_settings} '
            f'{optimization} -T {self._num_cores} {keepdirs}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _run_optimization(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`.Molecule`
            The optimized molecule.

        opt_complete : :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """

        xyz = 'input_structure_ff.xyz'
        out_file = 'crest_ff.output'
        mol.write(xyz)
        self._run_crest(xyz=xyz, out_file=out_file)

        # Check if the optimization is complete.
        output_xyzs = [
            'crest_best.xyz',
            'crest_conformers.xyz',
        ]
        opt_complete = self._is_complete(out_file, output_xyzs)
        mol = mol.with_structure_from_file(output_xyzs[0])

        return mol, opt_complete

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

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            mol, complete = self._run_optimization(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            logging.warning(f'CREST run is incomplete for {mol}.')

        return mol
