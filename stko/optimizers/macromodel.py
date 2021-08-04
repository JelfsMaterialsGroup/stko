"""
MacroModel Optimizers
=====================

Performs optimizations using MacroModel

"""


import logging
from uuid import uuid4

import rdkit.Chem.AllChem as rdkit

from ..packages.macromodel import (
    MacroModel,
    MacroModelInputError,
)
from ..utilities import (
    MAEExtractor,
    mol_from_mae_file,
    move_generated_macromodel_files,
)
from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class MacroModelOptimizer(Optimizer, MacroModel):
    """
    Base class for MacroModel optimizations.
    """

    def __init__(
        self,
        macromodel_path,
        output_dir,
        timeout,
        force_field,
        maximum_iterations,
        minimum_gradient,
    ):
        """
        Initialize a :class:`MacroModelOptimizer` object.

        Parameters
        ----------
        macromodel_path : :class:`str`
            The full path of the Schrodinger suite within the user's
            machine. For example, on a Linux machine this may be
            something like ``'/opt/schrodinger2017-2'``.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        timeout : :class:`float`, optional
            The amount in seconds the optimization is allowed to run
            before being terminated. ``None`` means there is no
            timeout.

        force_field : :class:`int`, optional
            The number of the force field to be used.
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
            |  16 | OPLS3/3e/4        |
            +------------+------------+

        minimum_gradient: : class: `float`
            The gradient at which optimization is stopped.
            Cannot be less than ``0.0001``.

        maximum_iterations: : class: `int`
            The maximum number of iterations done during the
            optimization. Cannot be more than ``999999``.
        """
        self._force_field = force_field
        self._minimum_gradient = minimum_gradient
        self._maximum_iterations = maximum_iterations
        super().__init__(
            macromodel_path=macromodel_path,
            output_dir=output_dir,
            timeout=timeout,
        )


class MacroModelForceField(MacroModelOptimizer):
    """
    Uses MacroModel force fields to optimize molecules.

    Examples
    --------
    Optimisation of any: class: `stk.Molecule` is possible with
    `restricted = False`.

    .. code-block: : python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCCN')
        optimizer = stko.MacroModelForceField(
            macromodel_path='/path/to/macromodel/',
        )
        mol = optimizer.optimize(mol)

    Optimisation of `long bonds` only within
    :class:`stk.ConstructedMolecule` is possible with
    restricted=True`.
    Generally, this means only bonds created during the construction
    process will be optimized,
    and those belonging to building blocks will be fixed.
    If the molecule is not a `ConstructedMolecule`, no positions will
    be optimized.

    .. code-block: : python

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
        optimizer = stko.MacroModelForceField(
            macromodel_path='/path/to/macromodel/',
            restricted=True,
        )
        polymer = optimizer.optimize(polymer)

    """

    def __init__(
        self,
        macromodel_path,
        output_dir=None,
        restricted=False,
        timeout=None,
        force_field=16,
        maximum_iterations=2500,
        minimum_gradient=0.05,
    ):
        """
        Initialize a: class: `MacroModelForceField` object.

        Parameters
        ----------
        macromodel_path: : class: `str`
            The full path of the Schrodinger suite within the user's
            machine. For example, on a Linux machine this may be
            something like ``'/opt/schrodinger2017-2'``.

        output_dir: : class: `str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            : func: `uuid.uuid4` is used.

        restricted: : class: `bool`, optional
            If ``True`` then an optimization is performed only on bonds
            not associated with building block IDs. These bonds may not
            correspond to the bonds formed between bonder atoms.
            If ``False`` then all bonds are optimized.

        timeout: : class: `float`, optional
            The amount in seconds the optimization is allowed to run
            before being terminated. ``None`` means there is no
            timeout.

        force_field: : class: `int`, optional
            The number of the force field to be used.
            Force field arguments can be the following:
            +------------+------------+
            |  1 | MM2                |
            +------------+------------+
            |  2 | MM3                |
            +------------+------------+
            |  3 | AMBER              |
            +------------+------------+
            |  4 | AMBER94            |
            +------------+------------+
            |  5 | OPLSA              |
            +------------+------------+
            |  10 | MMFF94 and MMFF94s|
            +------------+------------+
            |  14 | OPLS_2005         |
            +------------+------------+
            |  16 | OPLS3/3e/4        |
            +------------+------------+

        maximum_iterations: : class: `int`, optional
            The maximum number of iterations done during the
            optimization. Cannot be more than ``999999``.

        minimum_gradient: : class: `float`, optional
            The gradient at which optimization is stopped.
            Cannot be less than ``0.0001``.

        """
        self._check_params(
            minimum_gradient=minimum_gradient,
            maximum_iterations=maximum_iterations
        )
        self._restricted = restricted
        super().__init__(
            macromodel_path=macromodel_path,
            output_dir=output_dir,
            force_field=force_field,
            maximum_iterations=maximum_iterations,
            minimum_gradient=minimum_gradient,
            timeout=timeout,
        )

    @ staticmethod
    def _check_params(minimum_gradient, maximum_iterations):
        """
        Check if the optimization parameters are valid for MacroModel.

        Parameters
        ----------
        minimum_gradient: : class: `float`
            The gradient at which optimization is stopped.
            Cannot be less than ``0.0001``.

        maximum_iterations: : class: `int`
            The maximum number of iterations done during the
            optimization. Cannot be more than ``999999``.

        Returns
        -------
        None: : class: `NoneType`

        Raises
        ------
        : class: `.MacroModelInputError`
            If the parameters cannot be converted into a valid ``.com``
            file entry.

        """

        if minimum_gradient < 0.0001:
            raise MacroModelInputError(
                'Convergence gradient (< 0.0001) is too small.'
            )

        if maximum_iterations > 999999:
            raise MacroModelInputError(
                'Number of iterations (> 999999) is too high.'
            )

    def _generate_com(self, mol, run_name):
        """
        Create a ``.com`` file for a MacroModel optimization.

        The created ``.com`` file fixes all bond parameters which were
        not added by: meth: `~.Topology.construct`. This means all bond
        distances, bond angles and torsional angles are fixed, except
        for cases where it involves a bond added by
        : meth: `.Topology.construct`.

        This fixing is implemented by creating a ``.com`` file with
        various "FX" commands written within its body.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule which is to be optimized.

        run_name: : class: `str`
            The name of the run. The files generated by this run will
            have this name.

        Returns
        -------
        None: : class: `NoneType`

        """

        logger.debug(f'Creating .com file for "{mol}".')

        # This is the body of the ``.com`` file. The line that begins
        # and ends with exclamation lines is replaced with the various
        # commands that fix bond distances and angles.
        line1 = ('FFLD', self._force_field, 1, 0, 0, 1, 0, 0, 0)
        line2 = ('BGIN', 0, 0, 0, 0, 0, 0, 0, 0)
        line3 = ('READ', 0, 0, 0, 0, 0, 0, 0, 0)
        line4 = ('CONV', 2, 0, 0, 0, self._minimum_gradient, 0, 0, 0)
        line5 = ('MINI', 1, 0, self._maximum_iterations, 0, 0, 0, 0, 0)
        line6 = ('END', 0, 1, 0, 0, 0, 0, 0, 0)

        com_block = "\n".join([
            self._get_com_line(*line1),
            self._get_com_line(*line2),
            self._get_com_line(*line3),
            '!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!',
            self._get_com_line(*line4),
            self._get_com_line(*line5),
            self._get_com_line(*line6)
        ])

        # If `restricted` is ``False`` do not add a fix block.
        if not self._restricted:
            com_block = com_block.replace(
                "!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!\n",
                ''
            )
        else:
            # This function adds all the lines which fix bond distances
            # and angles into com_block.
            com_block = self._fix_params(mol, com_block)

        # Writes the .com file.
        with open(f'{run_name}.com', 'w') as com:
            # The first line holds the .mae file containing the
            # molecule to be optimized.
            com.write(f'{run_name}.mae\n')
            # The second line holds the name of the output file of the
            # optimization.
            com.write(f'{run_name}-out.maegz\n')
            # Next is the body of the .com file.
            com.write(com_block)

    def optimize(self, mol):
        """
        Optimize a molecule.

        Parameters
        ----------
        mol: : class: `stk.ConstructedMolecule`
            The molecule to be optimized.

        Returns
        -------
        mol: : class: `stk.ConstructedMolecule`
            The optimized molecule.

        """

        run_name = str(uuid4().int)
        if self._output_dir is None:
            output_dir = run_name
        else:
            output_dir = self._output_dir

        mol_path = f'{run_name}.mol'
        mae_path = f'{run_name}.mae'
        # First write a .mol file of the molecule.
        mol.write(mol_path)
        # MacroModel requires a ``.mae`` file as input.
        self._run_structconvert(mol_path, mae_path)
        # generate the ``.com`` file for the MacroModel run.
        self._generate_com(mol, run_name)
        # Run the optimization.
        self._run_bmin(mol, run_name)
        # Get the ``.maegz`` optimization output to a ``.mae``.
        self._convert_maegz_to_mae(run_name)
        rdkit_opt_mol = mol_from_mae_file(mae_path)
        mol = mol.with_position_matrix(
            rdkit_opt_mol.GetConformer().GetPositions()
        )
        move_generated_macromodel_files(run_name, output_dir)
        return mol

    def _fix_distances(self, mol, fix_block):
        """
        Add lines fixing bond distances to ``.com`` body.

        Only bond distances which do not involve bonds created during
        construction are fixed.

        Parameters
        ----------
        mol: : class: `stk.Molecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """
        # Identify bonds created by ``stk`` by checking if the
        # ``stk.BuildingBlock`` associated with the bond is ``None``.
        bonder_ids = set(
            atom_id
            for bond_info in mol.get_bond_infos()
            if bond_info.get_building_block() is None
            for atom_id in (
                bond_info.get_bond().get_atom1().get_id(),
                bond_info.get_bond().get_atom2().get_id()
            )
        )

        # Go through all the bonds in the ``stk.Molecule`` . If an
        # ``stk.BuildingBlock`` associated with the bond is ``None``,
        # add a fix line to the ``fix_block``.
        # If the bond is connected to an atom present in `bonder_ids`,
        # go to the next bond. This is because a bond between atoms,
        #  whose IDs are present in `bonder_ids`,
        # was added during construction and should therefore not
        # be fixed.
        for bond in mol.get_bonds():
            if (
                bond.get_atom1().get_id() in bonder_ids
                and
                bond.get_atom2().get_id() in bonder_ids
            ):
                continue

            atom1_id = bond.get_atom1().get_id()
            atom2_id = bond.get_atom2().get_id()
            # Make sure that the indices are increased by 1 in the .mae
            # file.
            atom1_id += 1
            atom2_id += 1
            args = ('FXDI', atom1_id, atom2_id, 0, 0, 99999, 0, 0, 0)
            fix_block += self._get_com_line(*args)
            fix_block += '\n'

        return fix_block

    def _fix_bond_angles(self, mol, fix_block):
        """
        Add lines fixing bond angles to the ``.com`` body.

        All bond angles of the molecule are fixed.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """

        paths = rdkit.FindAllPathsOfLengthN(
            mol=mol.to_rdkit_mol(),
            length=3,
            useBonds=False,
            useHs=True
        )
        for atom_ids in paths:
            atom_ids = [i+1 for i in atom_ids]
            args = ('FXBA', *atom_ids, 99999, 0, 0, 0, 0)
            fix_block += self._get_com_line(*args)
            fix_block += '\n'

        return fix_block

    def _fix_torsional_angles(self, mol, fix_block):
        """
        Add lines fixing torsional bond angles to the ``.com`` body.

        All torsional angles of the molecule are fixed.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """

        paths = rdkit.FindAllPathsOfLengthN(
            mol=mol.to_rdkit_mol(),
            length=4,
            useBonds=False,
            useHs=True
        )
        for atom_ids in paths:
            atom_ids = [i+1 for i in atom_ids]
            args = ('FXTA', *atom_ids, 99999, 361, 0, 0)
            fix_block += self._get_com_line(*args)
            fix_block += '\n'

        return fix_block


class MacroModelMD(MacroModelOptimizer):
    """
    Runs a molecular dynamics conformer search using MacroModel.

    Examples
    --------
    Molecular dynamics can be run on any: class: `stk.Molecule` using
    this class. Restrictions can be applied, but are not by default.
    This class collects a series of conformers from the trajectory,
    optimises them, then returns the lowest energy conformer.

    .. code-block: : python

        import stk
        import stko

        mol = stk.BuildingBlock('NCCCN')
        optimizer = stko.MacroModelMD(
            macromodel_path='/path/to/macromodel/',
            conformers=40,
        )
        mol = optimizer.optimize(mol)

    """

    def __init__(
        self,
        macromodel_path,
        output_dir=None,
        timeout=None,
        force_field=16,
        temperature=750,
        conformers=50,
        time_step=1,
        eq_time=10,
        simulation_time=200,
        maximum_iterations=2500,
        minimum_gradient=0.05,
        restricted_bonds=None,
        restricted_bond_angles=None,
        restricted_torsional_angles=None,
    ):
        """
        Initialize a: class: `.MacroModelMD` instance.

        Parameters
        ----------
        macromodel_path: : class: `str`
            The full path of the Schrodinger suite within the user's
            machine. For example, on a Linux machine this may be
            something like ``'/opt/schrodinger2017-2'``.

        output_dir: : class: `str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            : func: `uuid.uuid4` is used.

        timeout: : class: `float`, optional
            The amount in seconds the MD is allowed to run before
            being terminated. ``None`` means there is no timeout.

        force_field: : class: `int`, optional
            The number of the force field to be used.
            Force field arguments can be the following:
            +------------+------------+
            |  1 | MM2                |
            +------------+------------+
            |  2 | MM3                |
            +------------+------------+
            |  3 | AMBER              |
            +------------+------------+
            |  4 | AMBER94            |
            +------------+------------+
            |  5 | OPLSA              |
            +------------+------------+
            |  10 | MMFF94 and MMFF94s|
            +------------+------------+
            |  14 | OPLS_2005         |
            +------------+------------+
            |  16 | OPLS3/3e/4        |
            +------------+------------+


        temperature: : class: `float`, optional
            The temperature in Kelvin at which the MD is run.
            Cannot be more than ``99999.99``.

        conformers': : class: `int`, optional
            The number of conformers sampled and optimized from the MD.
            Cannot be more than ``9999``.

        simulation_time: : class: `float`, optional
            The simulation time in ``ps`` of the MD.
            Cannot be more than ``999999.99``.

        time_step: : class: `float`, optional
            The time step in ``fs`` for the MD.
            Cannot be more than ``99999.99``.

        eq_time: : class: `float`, optional
            The equilibration time in ``ps`` before the MD is run.
            Cannot be more than ``999999.99``.

        maximum_iterations: : class: `int`, optional
            The maximum number of iterations done during the
            optimization. Cannot be more than ``999999``.

        minimum_gradient: : class: `float`, optional
            The gradient at which optimization is stopped.
            Cannot be less than ``0.0001``.

        restricted_bonds: : class: `set`, optional
            A: class: `set` of the form

            .. code-block: : python

                restricted_bonds = {
                    frozenset((0, 10)),
                    frozenset((3, 14)),
                    frozenset((5, 6))
                }

            Where each: class: `frozenset` defines which bonds should
            have a fixed length via the atom ids of atoms in the bond.

        restricted_bond_angles: : class: `set`, optional
            A: class: `set` of the form

            .. code-block: : python

                restricted_bonds = {
                    frozenset((0, 10, 12)),
                    frozenset((3, 14, 7)),
                    frozenset((5, 8, 2))
                }

            Where each: class: `frozenset` defines which bond angles
            should have a fixed size via the atom ids of atoms in the
            bond angle.

        restricted_torsional_angles: : class: `set`, optional
            A: class: `set` of the form

            .. code-block: : python

                restricted_bonds = {
                    frozenset((0, 10, 12, 3)),
                    frozenset((3, 14, 7, 4)),
                    frozenset((5, 8, 2, 9))
                }

            Where each: class: `frozenset` defines which torsional
            angles should have a fixed size via the atom ids of atoms
            in the torsional angle.

        """

        if restricted_bonds is None:
            restricted_bonds = set()
        if restricted_bond_angles is None:
            restricted_bond_angles = set()
        if restricted_torsional_angles is None:
            restricted_torsional_angles = set()

        self._check_params(
            temperature=temperature,
            conformers=conformers,
            simulation_time=simulation_time,
            time_step=time_step,
            eq_time=eq_time,
            minimum_gradient=minimum_gradient,
            maximum_iterations=maximum_iterations
        )

        self._temperature = temperature
        self._conformers = conformers
        self._time_step = time_step
        self._eq_time = eq_time
        self._simulation_time = simulation_time
        self._restricted_bonds = restricted_bonds
        self._restricted_bond_angles = restricted_bond_angles
        self._restricted_torsional_angles = restricted_torsional_angles

        # Negative simulation time is interpreted as times 100 ps.
        if simulation_time > 99999.99:
            self._sim_time = -simulation_time/100
        else:
            self._sim_time = simulation_time

        # Negative equilibration time is interpreted as times 100 ps.
        if eq_time > 99999.99:
            self._eq_time = -eq_time/100
        else:
            self._eq_time = eq_time

        super().__init__(
            macromodel_path=macromodel_path,
            output_dir=output_dir,
            timeout=timeout,
            force_field=force_field,
            maximum_iterations=maximum_iterations,
            minimum_gradient=minimum_gradient,
        )

    @ staticmethod
    def _check_params(
        temperature,
        conformers,
        simulation_time,
        time_step,
        eq_time,
        minimum_gradient,
        maximum_iterations
    ):
        """
        Check if the optimization parameters are valid for MacroModel.

        Parameters
        ----------
        temperature: : class: `float`
            The temperature in Kelvin at which the MD is run.
            Cannot be more than ``99999.99``.

        conformers': : class: `int`
            The number of conformers sampled and optimized from the MD.
            Cannot be more than ``9999``.

        simulation_time: : class: `float`
            The simulation time in ``ps`` of the MD.
            Cannot be more than ``999999.99``.

        time_step: : class: `float`
            The time step in ``fs`` for the MD.
            Cannot be more than ``99999.99``.

        eq_time: : class: `float`
            The equilibriation time in ``ps`` before the MD is run.
            Cannot be more than ``999999.99``.

        minimum_gradient: : class: `float`
            The gradient at which optimization is stopped.
            Cannot be less than ``0.0001``.

        maximum_iterations: : class: `int`
            The maximum number of iterations done during the
            optimization. Cannot be more than ``999999``.

        Returns
        -------
        None: : class: `NoneType`

        Raises
        ------
        : class: `.MacroModelInputError`
            If the parameters cannot be converted into a valid ``.com``
            file entry.

        """

        if temperature > 99999.99:
            raise MacroModelInputError(
                'Supplied temperature (> 99999 K) is too high.'
            )

        if conformers > 9999:
            raise MacroModelInputError(
                'Supplied number of conformers (> 9999) is too high.'
            )

        if simulation_time > 999999.99:
            raise MacroModelInputError(
                'Supplied simulation time (> 999999 ps) is too long.'
            )

        if time_step > 99999.99:
            raise MacroModelInputError(
                'Supplied time step (> 99999 fs) is too high.'
            )

        if eq_time > 999999.99:
            raise MacroModelInputError(
                'Supplied eq time (> 999999 ps) is too long.'
            )

        if minimum_gradient < 0.0001:
            raise MacroModelInputError(
                'Convergence gradient (< 0.0001) is too small.'
            )

        if maximum_iterations > 999999:
            raise MacroModelInputError(
                'Number of iterations (> 999999) is too high.'
            )

    def _generate_com(self, mol, run_name):
        """
        Create a ``.com`` file for a MacroModel optimization.

        The created ``.com`` file fixes all bond parameters which were
        not added by: meth: `~.Topology.construct`. This means all bond
        distances, bond angles and torsional angles are fixed, except
        for cases where it involves a bond added by
        : meth: `~.Topology.construct`.

        This fixing is implemented by creating a ``.com`` file with
        various "FX" commands written within its body.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule which is to be optimized.

        Returns
        -------
        None: : class: `NoneType`

        """

        logger.debug(f'Creating .com file for "{mol}".')

        # Define some short aliases to keep the following lines neat.
        temp = self._temperature
        sim_time = self._sim_time
        tstep = self._time_step
        eq_time = self._eq_time

        line1 = ('FFLD', self._force_field, 1, 0, 0, 1, 0, 0, 0)
        line2 = ('READ', 0, 0, 0, 0, 0, 0, 0, 0)
        line3 = ('MDIT', 0, 0, 0, 0, temp, 0, 0, 0)
        line4 = ('MDYN', 0, 0, 0, 0, tstep, eq_time, temp, 0)
        line5 = ('MDSA', self._conformers, 0, 0, 0, 0, 0, 1, 0)
        line6 = ('MDYN', 1, 0, 0, 0, tstep, sim_time, temp, 0)
        line7 = ('WRIT', 0, 0, 0, 0, 0, 0, 0, 0)
        line8 = ('RWND', 0, 1, 0, 0, 0, 0, 0, 0)
        line9 = ('BGIN', 0, 0, 0, 0, 0, 0, 0, 0)
        line10 = ('READ', -2, 0, 0, 0, 0, 0, 0, 0)
        line11 = ('CONV', 2, 0, 0, 0, self._minimum_gradient, 0, 0, 0)
        line12 = (
            'MINI', 1, 0, self._maximum_iterations,
            0, 0, 0, 0, 0,
        )
        line13 = ('END', 0, 1, 0, 0, 0, 0, 0, 0)

        com_block = "\n".join([
            self._get_com_line(*line1),
            self._get_com_line(*line2),
            '!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!',
            self._get_com_line(*line3),
            self._get_com_line(*line4),
            self._get_com_line(*line5),
            self._get_com_line(*line6),
            self._get_com_line(*line7),
            self._get_com_line(*line8),
            self._get_com_line(*line9),
            self._get_com_line(*line10),
            '!!!BLOCK_OF_FIXED_PARAMETERS_COMES_HERE!!!',
            self._get_com_line(*line11),
            self._get_com_line(*line12),
            self._get_com_line(*line13),
        ])

        com_block = self._fix_params(mol, com_block)

        # Generate the com file containing the info for the run
        with open(f'{run_name}.com', 'w') as com:
            # name of the macromodel file
            com.write(f'{run_name}.mae\n')
            # name of the output file
            com.write(f'{run_name}-out.maegz\n')
            # details of the macromodel run
            com.write(com_block)

    def optimize(self, mol):
        """
        Optimize a molecule.

        Parameters
        ----------
        mol: : class: `.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol: : class: `.Molecule`
            The molecule to be optimized.

        """

        run_name = str(uuid4().int)
        if self._output_dir is None:
            output_dir = run_name
        else:
            output_dir = self._output_dir

        mol_path = f'{run_name}.mol'

        # First write a .mol file of the molecule.
        mol.write(mol_path)
        # MacroModel requires a ``.mae`` file as input.
        self._run_structconvert(mol_path, f'{run_name}.mae')
        # Generate the ``.com`` file for the MacroModel MD run.
        self._generate_com(mol, run_name)
        # Run the optimization.
        self._run_bmin(mol, run_name)
        # Extract the lowest energy conformer into its own .mae file.
        conformer_mae = MAEExtractor(run_name).path
        rdkit_opt_mol = mol_from_mae_file(conformer_mae)
        mol = mol.with_position_matrix(
            rdkit_opt_mol.GetConformer().GetPositions()
        )
        move_generated_macromodel_files(run_name, output_dir)
        return mol

    def _fix_distances(self, mol, fix_block):
        """
        Add lines fixing bond distances to ``.com`` body.

        Parameters
        ----------
        mol: : class: `stk.ConstructedMolecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """

        # Go through all the bonds in the rdkit molecule. If the bond
        # is not between bonder atoms add a fix line to the
        # ``fix_block``. If the bond does invovle two bonder atoms go
        # to the next bond. This is because a bond between 2 bonder
        # atoms was added during construction and should therefore not
        # be fixed.
        for bond in mol.get_bonds():
            bond_key = frozenset(
                (
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id()
                )
            )
            atom1_id = bond.get_atom1().get_id()
            atom2_id = bond.get_atom2().get_id()
            if (bond_key not in self._restricted_bonds):
                continue

            # Make sure that the indices are increased by 1 in the .mae
            # file.
            atom1_id += 1
            atom2_id += 1
            args = ('FXDI', atom1_id, atom2_id, 0, 0, 99999, 0, 0, 0)
            fix_block += self._get_com_line(*args)
            fix_block += '\n'

        return fix_block

    def _fix_bond_angles(self, mol, fix_block):
        """
        Add lines fixing bond angles to the ``.com`` body.

        Parameters
        ----------
        mol: : class: `stk.ConstructedMolecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """

        paths = rdkit.FindAllPathsOfLengthN(
            mol=mol.to_rdkit_mol(),
            length=3,
            useBonds=False,
            useHs=True
        )
        for atom_ids in paths:
            if frozenset(atom_ids) in self._restricted_bond_angles:
                atom_ids = [i+1 for i in atom_ids]
                args = ('FXBA', *atom_ids, 99999, 0, 0, 0, 0)
                fix_block += self._get_com_line(*args)
                fix_block += '\n'

        return fix_block

    def _fix_torsional_angles(self, mol, fix_block):
        """
        Add lines fixing torsional bond angles to the ``.com`` body.

        Parameters
        ----------
        mol: : class: `stk.ConstructedMolecule`
            The molecule to be optimized.

        fix_block: : class: `str`
            A string holding fix commands in the ``.com`` file.

        Returns
        -------
        : class: `str`
            A string holding fix commands in the ``.com`` file.

        """

        paths = rdkit.FindAllPathsOfLengthN(
            mol=mol.to_rdkit_mol(),
            length=4,
            useBonds=False,
            useHs=True
        )
        # Apply the fix.
        for atom_ids in paths:
            if frozenset(atom_ids) in (
                self._restricted_torsional_angles
            ):
                atom_ids = [i+1 for i in atom_ids]
                args = ('FXTA', *atom_ids, 99999, 361, 0, 0)
                fix_block += self._get_com_line(*args)
                fix_block += '\n'

        return fix_block
