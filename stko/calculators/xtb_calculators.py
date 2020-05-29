"""
XTB Calculators
===============

#. :class:`.XTBEnergy`

Wrappers for calculators within the :mod:`xtb` code.

"""

import logging
import os
import shutil
import uuid
import subprocess as sp

from ..base_calculator import Calculator
from ..utilities import (
    is_valid_xtb_solvent,
    XTBInvalidSolventError,
    XTBExtractor
)

logger = logging.getLogger(__name__)


class XTBEnergy(Calculator):
    """
    Uses GFN-xTB [1]_ to calculate energy and other properties.

    By default, :meth:`get_energy` will extract other properties of the
    :class:`.Molecule` passed to :meth:`get_energy`, which
    will be saved in the attributes of :class:`.XTBEnergy`.

    Notes
    -----
    When running :meth:`get_energy`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTBEnergy` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Attributes
    ----------
    total_energies : :class:`dict`
        :class:`dict` of the total energy of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    homo_lumo_gaps : :class:`dict`
        :class:`dict` of the HOMO-LUMO gap of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    fermi_levels : :class:`dict`
        :class:`dict` of the Fermi level of each :class:`.Molecule`
        passed to :meth:`get_energy`.

    homo_lumo_orbitals : :class:`dict`
        :class:`dict` of the HOMO-LUMO orbital properties of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qonly_dipole_moments : :class:`dict`
        :class:`dict` of the q only dipole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    full_dipole_moments : :class:`dict`
        :class:`dict` of the full dipole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qonly_quadrupole_moments : :class:`dict`
        :class:`dict` of the q only quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    qdip_quadrupole_moments : :class:`dict`
        :class:`dict` of the q+dip quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    full_quadrupole_moments : :class:`dict`
        :class:`dict` of the full quadrupole moment of
        each :class:`.Molecule` passed to :meth:`get_energy`.

    total_free_energies : :class:`dict`
        :class:`dict` of the total free energy of
        each :class:`.Molecule` passed to :meth:`get_energy`.
        This is empty if :attr:`calculate_free_energy` is ``False``.

    frequencies : :class:`dict`
        :class:`dict` of the vibrational frequencies of
        each :class:`.Molecule` passed to :meth:`get_energy`.
        This is empty if `calculate_free_energy` was ``False``.

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
        opt = stko.OptimizerSequence(
            stko.UFF(),
            stko.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
        )
        polymer = opt.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = stk.XTBEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True
        )

        p_total_energy = xtb.get_energy(polymer)

        # Extract properties from the energy calculator for a given
        # molecule.
        homo_lumo_gap = xtb.homo_lumo_gaps[polymer]
        fermi_levels = xtb.fermi_levels[polymer]
        homo_lumo_orbitals = xtb.homo_lumo_orbitals[polymer]
        full_dipole_moments = xtb.full_dipole_moments[polymer]
        full_quadrupole_moments = xtb.full_quadrupole_moments[polymer]

        # The total energy can be extracted at any point from the
        # calculator.
        total_energy = xtb.total_energies[polymer]

    If `calculate_free_energy` is ``True``, xTB performs a
    numerical Hessian calculation and calculates the total free energy
    and vibrational frequencies of a molecule. It is recommended that a
    well optimized structure be used as input for these calculations

    .. code-block:: python

        # Optimize the constructed molecule so that it has a
        # reasonable structure.
        optimizer = stk.OptimizerSequence(
            stk.ETKDG(),
            stk.XTB(
                xtb_path='/opt/gfnxtb/xtb',
                unlimited_memory=True,
                opt_level='verytight'
            )
        )
        polymer = optimizer.optimize(polymer)

        # Calculate energy using GFN-xTB.
        xtb = stk.XTBEnergy(
            xtb_path='/opt/gfnxtb/xtb',
            unlimited_memory=True,
            calculate_free_energy=True
        )

        p_total_energy = xtb.get_energy(polymer)

        # Extract properties from the energy calculator for a given
        # molecule.
        p_total_free_energy = xtb.total_energies[polymer]
        p_frequencies = xtb.frequencies[polymer]

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html

    """
    def __init__(
        self,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        num_cores=1,
        calculate_free_energy=False,
        electronic_temperature=300,
        solvent=None,
        solvent_grid='normal',
        charge=0,
        num_unpaired_electrons=0,
        unlimited_memory=False,
    ):
        """
        Initializes a :class:`XTBEnergy` instance.

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
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        calculate_free_energy : :class:`bool`, optional
            Whether to calculate the total free energy and vibrational
            frequencies. Setting this to ``True`` can drastically
            increase calculation time and memory requirements.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

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
            If ``True`` :meth:`energy` will be run without constraints
            on the stack size. If memory issues are encountered, this
            should be ``True``, however this may raise issues on
            clusters.

        """
        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    f'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(gfn_version, solvent):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} is invalid for ',
                    f'version {gfn_version!r}.'
                )

        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
        self._num_cores = str(num_cores)
        self._calculate_free_energy = calculate_free_energy
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory

        self.total_energies = {}
        self.homo_lumo_gaps = {}
        self.fermi_levels = {}
        self.homo_lumo_orbitals = {}
        self.qonly_dipole_moments = {}
        self.full_dipole_moments = {}
        self.qonly_quadrupole_moments = {}
        self.qdip_quadrupole_moments = {}
        self.full_quadrupole_moments = {}
        self.total_free_energies = {}
        self.frequencies = {}

    def _get_properties(self, mol, output_file):
        """
        Extracts properties from a GFN-xTB energy calculation.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The :class:`.Molecule` whose energy was calculated.

        output_file : :class: `str`
            Name of the output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Get properties from output_file.
        xtbext = XTBExtractor(output_file=output_file)

        self.total_energies[mol] = xtbext.total_energy
        self.homo_lumo_gaps[mol] = xtbext.homo_lumo_gap
        self.fermi_levels[mol] = xtbext.fermi_level
        self.homo_lumo_orbitals[mol] = xtbext.homo_lumo_occ
        self.qonly_dipole_moments[mol] = xtbext.qonly_dipole_moment
        self.full_dipole_moments[mol] = xtbext.full_dipole_moment
        self.qonly_quadrupole_moments[mol] = \
            xtbext.qonly_quadrupole_moment
        self.qdip_quadrupole_moments[mol] = \
            xtbext.qdip_quadrupole_moment
        self.full_quadrupole_moments[mol] = \
            xtbext.full_quadrupole_moment
        if self._calculate_free_energy:
            self.total_free_energies[mol] = xtbext.total_free_energy
            self.frequencies[mol] = xtbext.frequencies

    def _run_xtb(self, xyz, out_file):
        """
        Runs GFN-xTB.

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
            solvent = f'--gbsa {self._solvent} {self._solvent_grid}'
        else:
            solvent = ''

        if self._calculate_free_energy:
            calc_type = '--hess'
        else:
            calc_type = ''

        cmd = (
            f'{memory} {self._xtb_path} '
            f'{xyz} --gfn {self._gfn_version} '
            f'{calc_type} --parallel {self._num_cores} '
            f'--etemp {self._electronic_temperature} '
            f'{solvent} --chrg {self._charge} '
            f'--uhf {self._num_unpaired_electrons}'
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

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        xyz = os.path.join(output_dir, 'input_structure.xyz')
        out_file = os.path.join(output_dir, 'energy.output')
        mol.write(xyz)

        try:
            os.chdir(output_dir)
            self._run_xtb(xyz=xyz, out_file=out_file)
        finally:
            os.chdir(init_dir)

        self._get_properties(
            mol=mol,
            output_file=out_file
        )
        return self.total_energies[mol]
