import logging
import math
import os
import re
import shutil
import subprocess as sp
import uuid

import stk
from rdkit.Chem import AllChem as rdkit
from stko._internal.molecular.periodic.unitcell import UnitCell
from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.optimizers.utilities import (
    get_metal_atoms,
    get_metal_bonds,
    has_h_atom,
    has_metal_atom,
    to_rdkit_mol_without_metals,
)
from stko._internal.utilities.exceptions import (
    ExpectedMetalError,
    ForceFieldSetupError,
    OptimizerError,
    PathError,
)

logger = logging.getLogger(__name__)


class GulpUFFOptimizer(Optimizer):
    """Applies forcefield optimizers that can handle metal centres.

    Notes:
    -----
        By default, :meth:`optimize` will run an optimisation using the
        UFF4MOF. This forcefield requires some explicit metal atom
        definitions, which are determined by the user.

        This code was originally written for use with Gulp 5.1 on Linux and has
        not been officially tested on other versions and operating systems.
        Make sure to sanity check the output.

    Examples:
    --------
        While metal atoms are not required, UFF4MOF is useful because it
        encompasses almost all chemical environments commonly found in
        metal-organic structures. Better forcefields exist for purely
        organic molecules! An interface with GULP is provided, which takes
        the forcefield types assigned by RDKit for non-metal atoms and
        user defined forcefield types for metal atoms to perform geometry
        optimisations.

        .. code-block:: python

            import stk
            import stko
            from rdkit.Chem import AllChem as rdkit

            # Produce a Pd+2 atom with 4 functional groups.
            atom = rdkit.MolFromSmiles('[Pd+2]')
            atom.AddConformer(rdkit.Conformer(atom.GetNumAtoms()))
            palladium_atom = stk.BuildingBlock.init_from_rdkit_mol(atom)
            atom_0, = palladium_atom.get_atoms(0)
            palladium_atom = palladium_atom.with_functional_groups(
                (stk.SingleAtom(atom_0) for i in range(4))
            )

            # Build a building block with two functional groups using
            # the SmartsFunctionalGroupFactory.
            bb1 = stk.BuildingBlock(
                smiles=('C1=CC(=CC(=C1)C2=CN=CC=C2)C3=CN=CC=C3'),
                functional_groups=[
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[#6]~[#7X2]~[#6]',
                        bonders=(1, ),
                        deleters=(),
                    ),
                ],
            )

            # Build a metal-organic cage with dative bonds between
            # GenericFunctionalGroup and SingleAtom functional groups.
            cage = stk.ConstructedMolecule(
                stk.cage.M2L4Lantern(
                    building_blocks={
                        palladium_atom: (0, 1),
                        bb1: (2, 3, 4, 5)
                    },
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    )
                )
            )

            # Perform Gulp optimisation with UFF4MOF.
            # Use conjugate gradient method for a slower, but more stable
            # optimisation.

            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path='path/to/gulp',
                metal_FF={46: 'Pd4+2'},
                conjugate_gradient=True
            )

            # Assign the force field.
            gulp_opt.assign_FF(cage)
            # Run optimization.
            cage = gulp_opt.optimize(mol=cage)

    """

    def __init__(
        self,
        gulp_path: str,
        maxcyc: int = 1000,
        metal_FF: dict | None = None,
        metal_ligand_bond_order: str | None = None,
        conjugate_gradient: bool = False,
        output_dir: str | None = None,
    ):
        """Parameters
        gulp_path:
            Path to GULP executable.

        maxcyc:
            Set the maximum number of optimisation steps to use.
            Default in Gulp is 1000.

        metal_FF:
            Dictionary with metal atom forcefield assignments.
            Key: :class:`int` : atomic number.
            Value: :class:`str` : UFF4MOF forcefield type.

        metal_ligand_bond_order:
            Bond order to use for metal-ligand bonds. Defaults to
            `half`, but using `resonant` can increase the force
            constant for stronger metal-ligand interactions.

        conjugate_gradient:
            ``True`` to use Conjugate Graditent method.
            Defaults to ``False``

        output_dir:
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        """
        self._check_path(gulp_path)
        self._gulp_path = gulp_path
        self._maxcyc = maxcyc
        self._metal_FF = metal_FF
        self._metal_ligand_bond_order = (
            "half"
            if metal_ligand_bond_order is None
            else metal_ligand_bond_order
        )
        self._conjugate_gradient = conjugate_gradient
        self._output_dir = output_dir

    def _check_path(self, path: str) -> None:
        if not os.path.exists(path):
            raise PathError(f"GULP not found at {path}")

    def _add_atom_charge_flags(self, atom: rdkit.Atom, atomkey: str) -> str:
        """Add atom charge flags for forcefield.

        Code inspired by:
        https://github.com/rdkit/rdkit
        >   Code/GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        total_valence = rdkit.Atom.GetTotalValence(atom)
        atnum = int(atom.GetAtomicNum())

        # Go through element cases.
        # Mg.
        if atnum == 12:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Al.
        elif atnum == 13:
            if total_valence != 3:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Si.
        elif atnum == 14:
            if total_valence != 4:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # P.
        elif atnum == 15:
            if total_valence == 3:
                atomkey += "+3"
            elif total_valence == 5:
                atomkey += "+5"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # S.
        elif atnum == 16:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if hybrid != rdkit.HybridizationType.SP2:
                if total_valence == 2:
                    atomkey += "+2"
                elif total_valence == 4:
                    atomkey += "+4"
                elif total_valence == 6:
                    atomkey += "+6"
                else:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized charge state for "
                        f"atom: {atom.GetIdx}"
                    )
        # Zn.
        elif atnum == 30:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Ga.
        elif atnum == 31 or atnum == 33:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Se.
        elif atnum == 34 or atnum == 48:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # In.
        elif atnum == 49 or atnum == 51:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Te.
        elif atnum == 52 or atnum == 80:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Tl.
        elif atnum == 81 or atnum == 82 or atnum == 83:
            if total_valence == 3:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Po.
        elif atnum == 84:
            if total_valence == 2:
                atomkey += "+2"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Lanthanides.
        elif atnum >= 57 and atnum <= 71:
            if total_valence == 6:
                atomkey += "+3"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        return atomkey

    def _get_atom_label(self, atom: rdkit.Atom) -> str:
        """Get FF atom label.

        Code inspired by:
        https://github.com/rdkit/rdkit
        >   Code/GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        atnum = int(atom.GetAtomicNum())
        atomkey = atom.GetSymbol()
        if len(atomkey) == 1:
            atomkey += "_"

        table = rdkit.GetPeriodicTable()

        chk1 = rdkit.PeriodicTable.GetDefaultValence(table, atnum) == -1
        chk2 = rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 1
        chk3 = rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 7
        chk4 = chk2 and chk3
        if chk1 or chk4:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if atnum == 84:
                atomkey += "3"
                if hybrid != rdkit.HybridizationType.SP3:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            elif atnum == 80:
                atomkey += "1"
                if hybrid != rdkit.HybridizationType.SP:
                    raise ForceFieldSetupError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            elif hybrid == rdkit.HybridizationType.SP:
                atomkey += "1"
            elif hybrid == rdkit.HybridizationType.SP2:
                chk1a = rdkit.Atom.GetIsAromatic(atom)
                bonds = rdkit.Atom.GetBonds(atom)
                conjugated = False
                for bond in bonds:
                    if rdkit.Bond.GetIsConjugated(bond):
                        conjugated = True
                        break
                chk2a = conjugated
                chk3a = atnum in [6, 7, 8, 16]
                chk4a = chk1a or chk2a
                if chk4a and chk3a:
                    atomkey += "R"
                else:
                    atomkey += "2"
            elif hybrid == rdkit.HybridizationType.SP3:
                atomkey += "3"
            elif hybrid == rdkit.HybridizationType.SP3D:
                atomkey += "5"
            elif hybrid == rdkit.HybridizationType.SP3D2:
                atomkey += "6"
            else:
                raise ForceFieldSetupError(
                    f"UFFTYPER: Unrecognized hybridization for"
                    f" atom: {atom.GetIdx}"
                )
        atomkey = self._add_atom_charge_flags(atom, atomkey)
        return atomkey

    def _type_translator(self) -> dict[str, str]:
        type_translator: dict[str, str] = {}
        types = sorted(
            set(
                [  # type: ignore[type-var]
                    self.atom_labels[i][0] for i in self.atom_labels
                ]
            )
        )
        for t in types:
            if not t[1].isalpha():  # type:ignore[index]
                symb = t[0]  # type:ignore[index]
            else:
                symb = t[0:2]  # type:ignore[index]
            for i in range(1, 100):
                name = f"{symb}{i}"
                if name in type_translator.values():
                    continue
                else:
                    type_translator[t] = name  # type: ignore[index]
                    break

        return type_translator

    def _cell_section(self, unit_cell: UnitCell) -> str:
        cell_section = (
            "\ncell\n"
            f"{round(unit_cell.get_a(), 6)} "
            f"{round(unit_cell.get_b(), 6)} "
            f"{round(unit_cell.get_c(), 6)} "
            f"{round(unit_cell.get_alpha(), 6)} "
            f"{round(unit_cell.get_beta(), 6)} "
            f"{round(unit_cell.get_gamma(), 6)} "
            # No fixes.
            "0 0 0 0 0 0\n"
        )

        return cell_section

    def _position_section(
        self, mol: stk.Molecule, type_translator: dict
    ) -> str:
        position_section = "\ncartesian\n"
        for atom in mol.get_atoms():
            atom_type = type_translator[self.atom_labels[atom.get_id()][0]]
            position = mol.get_centroid(atom_ids=atom.get_id())
            posi_string = (
                f"{atom_type} core {round(position[0], 5)} "
                f"{round(position[1], 5)} {round(position[2], 5)}\n"
            )
            position_section += posi_string

        return position_section

    def _bond_section(
        self,
        mol: stk.Molecule,
        metal_atoms: list[stk.Atom],
    ) -> str:
        bond_section = "\n"
        for bond in mol.get_bonds():
            atom_types = [
                self.atom_labels[i.get_id()][0]
                for i in [bond.get_atom1(), bond.get_atom2()]
            ]

            # Set bond orders.
            if has_h_atom(bond):
                # H has bond order of 1.
                bond_type = ""
            elif has_metal_atom(bond, metal_atoms):
                bond_type = self._metal_ligand_bond_order
            elif (
                "_R" in atom_types[0]  # type:ignore[operator]
                and "_R" in atom_types[1]  # type:ignore[operator]
            ):
                bond_type = "resonant"
            elif bond.get_order() == 1:
                bond_type = ""
            elif bond.get_order() == 2:
                bond_type = "double"
            elif bond.get_order() == 3:
                bond_type = "triple"

            string = (
                f"connect {bond.get_atom1().get_id()+1} "
                f"{bond.get_atom2().get_id()+1} {bond_type}"
            )
            bond_section += string + "\n"

        return bond_section

    def _species_section(self, type_translator: dict) -> str:
        species_section = "\nspecies\n"
        for spec in type_translator:
            name = type_translator[spec]
            species_section += f"{name} {spec}\n"

        return species_section

    def _write_gulp_file(
        self,
        mol: stk.Molecule,
        metal_atoms: list[stk.Atom],
        in_file: str,
        output_xyz: str,
        unit_cell: UnitCell | None = None,
    ) -> None:
        type_translator = self._type_translator()

        top_line = "opti "

        if self._conjugate_gradient:
            top_line += "conj unit "

        if unit_cell is not None:
            # Constant pressure.
            top_line += "conp "
            cell_section = self._cell_section(unit_cell)
            # Output CIF.
            output_cif = output_xyz.replace("xyz", "cif")
            periodic_output = f"output cif {output_cif}\n"
        else:
            # Constant volume.
            top_line += "conv "
            cell_section = ""
            periodic_output = ""

        top_line += "noautobond fix molmec cartesian\n"

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = "\nlibrary uff4mof.lib\n"

        output_section = (
            "\n"
            f"maxcyc {self._maxcyc}\n"
            "terse inout potentials\n"
            "terse in cell\n"
            "terse in structure\n"
            "terse inout derivatives\n"
            f"output xyz {output_xyz}\n"
            f"{periodic_output}"
            # 'output movie xyz steps_.xyz\n'
        )

        with open(in_file, "w") as f:
            f.write(top_line)
            f.write(cell_section)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(output_section)

    def assign_FF(self, mol: stk.Molecule) -> None:
        """Assign forcefield types to molecule.

        Parameters
        ----------
            mol:
                The molecule to be optimized.

        """
        FutureWarning(
            "We have found some minor discrepancies in this "
            "assignment algorithm, which is based off rdkit code. "
            "Changes should come soon. This UFF optimisation should "
            " not be your final step! Due to this, some tests in "
            "test_uff_assign_ff.py have been muted."
        )

        metal_atoms = get_metal_atoms(mol)
        metal_ids = [i.get_id() for i in metal_atoms]

        if len(metal_ids) > 1 and self._metal_FF is None:
            raise ExpectedMetalError(
                "No metal FF provivded, but metal atoms were found ("
                f"{metal_atoms})"
            )

        metal_bonds, _ = get_metal_bonds(mol, metal_atoms)
        edit_mol = to_rdkit_mol_without_metals(
            mol=mol, metal_atoms=metal_atoms, metal_bonds=metal_bonds
        )

        # Get forcefield parameters.
        rdkit.SanitizeMol(edit_mol)
        self.atom_labels = {}

        for i in range(edit_mol.GetNumAtoms()):
            if i in metal_ids:
                self.atom_labels[i] = [None, "metal", None]
            else:
                atom = edit_mol.GetAtomWithIdx(i)
                atom_label = self._get_atom_label(atom)
                self.atom_labels[i] = [atom_label, None, None]

        # Write UFF4MOF specific forcefield parameters.
        # Metals.
        for atomid in self.atom_labels:
            if self.atom_labels[atomid][1] == "metal":
                (atom,) = mol.get_atoms(atomid)
                atom_no = atom.get_atomic_number()
                self.atom_labels[atomid][0] = self._metal_FF[  # type:ignore[index]
                    atom_no
                ]

    def _run_gulp(self, in_file: str, out_file: str) -> None:
        cmd = f"{self._gulp_path} < {in_file}"
        with open(out_file, "w") as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True,
            )

    def extract_final_energy(self, out_file: str) -> float:
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        with open(out_file) as f:
            for line in f.readlines():
                if "Final energy =" in line:
                    string = nums.search(line.rstrip())
                    return float(string.group(0))  # type: ignore[union-attr]

        raise OptimizerError(
            f'"Final energy =" not found in {out_file}, implying unsuccesful'
            " optimisation"
        )

    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
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

        in_file = "gulp_opt.gin"
        out_file = "gulp_opt.ginout"
        output_xyz = "gulp_opt.xyz"

        metal_atoms = get_metal_atoms(mol)

        try:
            # Write GULP file.
            self._write_gulp_file(
                mol=mol,
                metal_atoms=metal_atoms,
                in_file=in_file,
                output_xyz=output_xyz,
                unit_cell=None,
            )
            # Run.
            self._run_gulp(in_file, out_file)

            # Update from output.
            mol = mol.with_structure_from_file(output_xyz)

        finally:
            os.chdir(init_dir)

        return mol

    def p_optimize(
        self,
        mol: stk.Molecule,
        unit_cell: UnitCell,
    ) -> tuple[stk.Molecule, UnitCell]:
        """Optimize `mol` and `unit_cell`.

        Parameters
        ----------
            mol:
                The molecule to be optimized.

            unit_cell:
                The unit_cell to be optimized if optimization is periodic.

        Returns:
        -------
            mol:
                The optimized molecule.

            unit_cell:
                The optimized cell.

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

        in_file = "gulp_opt.gin"
        out_file = "gulp_opt.ginout"
        output_xyz = "gulp_opt.xyz"
        output_cif = output_xyz.replace("xyz", "cif")

        metal_atoms = get_metal_atoms(mol)

        try:
            # Write GULP file.
            self._write_gulp_file(
                mol=mol,
                metal_atoms=metal_atoms,
                in_file=in_file,
                output_xyz=output_xyz,
                unit_cell=unit_cell,
            )

            # Run.
            self._run_gulp(in_file, out_file)

            # Update from output.
            mol = mol.with_structure_from_file(output_xyz)
            unit_cell = unit_cell.with_cell_from_cif(output_cif)

        finally:
            os.chdir(init_dir)

        return mol, unit_cell


class GulpUFFMDOptimizer(GulpUFFOptimizer):
    """Applies forcefield MD that can handle metal centres.

    Notes:
    -----
        By default, :meth:`optimize` will run a MD run using the UFF4MOF.
        This forcefield requires some explicit metal atom definitions,
        which are determined by the user.

        This code was originally written for use with Gulp 5.1 on Linux
        and has not been officially tested on other versions and operating
        systems. Make sure to sanity check the output.


    Examples:
    --------
        Conformer searching is often useful, so we have provided an interface
        to MD simulations using GULP and UFF4MOF. A conformer search can be
        run at high temperature, where N conformers are extracted at constant
        intervals throughtout the simulation and optimized using UFF4MOF. The
        lowest energy conformer is returned. After these MD steps, it is
        crucial to reoptimize the resultant structure using a better
        forcefield or a more robust method!

        .. code-block:: python

            gulp_MD = stko.GulpUFFMDOptimizer(
                gulp_path='path/to/gulp',
                metal_FF={46: 'Pd4+2'},
                temperature=300,
                N_conformers=10,
                opt_conformers=True,
            )
            gulp_MD.assign_FF(cage)
            cage = gulp_MD.optimize(cage)

    """

    def __init__(
        self,
        gulp_path: str,
        metal_FF: dict[int, str] | None = None,
        metal_ligand_bond_order: str | None = None,
        output_dir: str | None = None,
        integrator: str = "stochastic",
        ensemble: str = "nvt",
        temperature: float = 300,
        equilbration: float = 1.0,
        production: float = 10.0,
        timestep: float = 1.0,
        N_conformers: int = 10,
        opt_conformers: bool = True,
        save_conformers: bool = False,
    ) -> None:
        """Parameters
        gulp_path:
            Path to GULP executable.

        metal_FF:
            Dictionary with metal atom forcefield assignments.
            Key: :class:`int` : atomic number.
            Value: :class:`str` : UFF4MOF forcefield type.

        metal_ligand_bond_order:
            Bond order to use for metal-ligand bonds. Defaults to
            `half`, but using `resonant` can increase the force
            constant for stronger metal-ligand interactions.

        output_dir:
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        integrator:
            Integrator for GULP to use.
            Defaults to 'stochastic'.

        ensemble:
            Ensemble for GULP to use.
            Defaults to 'nvt'.

        temperature:
            Temperature to run simulation at in Kelvin.
            Defaults to 300.

        equilbration:
            Time spent equilibrating system in ps.
            Defaults to 1.0.

        production:
            Time spent running production simulation of system in ps.
            Defaults to 10.0.

        timestep:
            Timestep of simulation in fs.
            Defaults to 1.0.

        N_conformers:
            Number of conformers to sample.
            Defaults to 10.

        opt_conformers:
            Whether or not to optimise each conformer using UFF4MOF.
            Defaults to ``True``.

        save_conformers:
            Whether or not to save to file each conformer.
            Defaults to ``False``.

        """
        self._check_path(gulp_path)
        self._gulp_path = gulp_path
        self._metal_FF = metal_FF
        self._metal_ligand_bond_order = (
            "half"
            if metal_ligand_bond_order is None
            else metal_ligand_bond_order
        )
        self._output_dir = output_dir
        self._integrator = integrator
        self._ensemble = ensemble
        self._temperature = temperature
        self._equilbration = float(equilbration)
        self._production = float(production)
        self._timestep = timestep
        self._N_conformers = N_conformers
        samples = float(self._production) / float(self._N_conformers)
        samples = math.floor(samples * 1000) / 1000
        self._sample = samples
        self._write = samples
        self._opt_conformers = opt_conformers
        self._save_conformers = save_conformers

    def _write_gulp_file(
        self,
        mol: stk.Molecule,
        metal_atoms: list[stk.Atom],
        in_file: str,
        output_traj: str,
        unit_cell: UnitCell | None = None,
    ) -> None:
        type_translator = self._type_translator()

        top_line = "md "

        if unit_cell is not None:
            # Constant pressure.
            top_line += "conp "
            cell_section = self._cell_section(unit_cell)
            # Output CIF.
            output_cif = output_traj.replace("xyz", "cif")
            periodic_output = f"output cif {output_cif}\n"
        else:
            # Constant volume.
            top_line += "conv "
            cell_section = ""
            periodic_output = ""

        top_line += "noautobond fix molmec cartesian\n"

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = "\nlibrary uff4mof.lib\n"

        md_section = (
            "\n"
            "mdmaxtemp 100000000\n"
            f"integrator {self._integrator}\n"
            f"ensemble {self._ensemble}\n"
            f"temperature {self._temperature}\n"
            f"equilibration {self._equilbration} ps\n"
            f"production {self._production} ps\n"
            f"timestep {self._timestep} fs\n"
            f"sample {self._sample} ps\n"
            f"write {self._write} ps\n"
            "\n"
        )

        output_section = (
            "\n"
            "terse inout potentials\n"
            "terse inout cell\n"
            "terse in structure\n"
            "terse inout derivatives\n"
            f"output trajectory ascii {output_traj}\n"
            f"{periodic_output}"
        )

        with open(in_file, "w") as f:
            f.write(top_line)
            f.write(cell_section)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(md_section)
            f.write(output_section)

    def _convert_traj_to_xyz(
        self,
        output_xyz: str,
        output_traj: str,
    ) -> tuple[list[str], dict[int, dict], list[str]]:
        # Get atom types from an existing xyz file.
        atom_types = []
        with open(output_xyz) as f:
            for line in f.readlines()[2:]:
                atom_types.append(line.rstrip().split(" ")[0])

        # Read in lines from trajectory file.
        with open(output_traj) as f:
            lines = f.readlines()

        # Split file using strings.
        timesteps = "".join(lines).split("#  Time/KE/E/T")[1:]
        trajectory_data = {}
        xyz_traj_lines = []
        for ts, cont in enumerate(timesteps):
            ts_data = {}
            time_section = cont.split("#  Coordinates\n")[0]
            coords_section = cont.split("#  Coordinates\n")[1].split(
                "#  Velocities\n"
            )[0]
            vels_section = cont.split("#  Velocities\n")[1].split(
                "#  Derivatives \n"
            )[0]
            derivs_section = cont.split("#  Derivatives \n")[1].split(
                "#  Site energies \n"
            )[0]
            sites_section = cont.split("#  Site energies \n")[1]

            ts_data["time"] = float(
                [i for i in time_section.strip().split(" ") if i][0]
            )
            ts_data["KE"] = float(
                [i for i in time_section.strip().split(" ") if i][1]
            )
            ts_data["E"] = float(
                [i for i in time_section.strip().split(" ") if i][2]
            )
            ts_data["T"] = float(
                [i for i in time_section.strip().split(" ") if i][3]
            )
            ts_data["coords"] = [  # type:ignore[assignment]
                [i for i in li.split(" ") if i]
                for li in coords_section.split("\n")[:-1]
            ]
            ts_data["vels"] = [  # type:ignore[assignment]
                [i for i in li.split(" ") if i]
                for li in vels_section.split("\n")[:-1]
            ]
            ts_data["derivs"] = [  # type:ignore[assignment]
                [i for i in li.split(" ") if i]
                for li in derivs_section.split("\n")[:-1]
            ]
            ts_data["sites"] = [  # type:ignore[assignment]
                [i for i in li.split(" ") if i]
                for li in sites_section.split("\n")[:-1]
            ]

            trajectory_data[ts] = ts_data

            # Write XYZ string for XYZ traj file.
            xyz_string = (
                f"{len(ts_data['coords'])}\n"  # type:ignore[arg-type]
                f"{ts_data['time']},{ts_data['KE']},"
                f"{ts_data['E']},{ts_data['T']}\n"
            )

            for i, coord in enumerate(  # type: ignore[var-annotated]
                ts_data["coords"]  # type: ignore[arg-type]
            ):
                site_E = ts_data["sites"][i][0]  # type: ignore[index]
                xyz_string += (
                    f"{atom_types[i]} {round(float(coord[0]), 5)} "
                    f"{round(float(coord[1]), 5)} "
                    f"{round(float(coord[2]), 5)} {site_E}\n"
                )

            xyz_traj_lines.append(xyz_string)

        return atom_types, trajectory_data, xyz_traj_lines

    def _write_conformer_xyz_file(
        self,
        ts: float,
        ts_data: dict,
        filename: str,
        atom_types: list[str],
    ) -> None:
        coords = ts_data["coords"]
        xyz_string = (
            f"{len(ts_data['coords'])}\n"
            f"{ts},{ts_data['time']},{ts_data['KE']},"
            f"{ts_data['E']},{ts_data['T']}\n"
        )

        for i, coord in enumerate(coords):
            xyz_string += (
                f"{atom_types[i]} {round(float(coord[0]), 5)} "
                f"{round(float(coord[1]), 5)} "
                f"{round(float(coord[2]), 5)}\n"
            )
        with open(filename, "w") as f:
            f.write(xyz_string)

    def _optimise_all_conformers(
        self,
        mol: stk.Molecule,
        trajectory_data: dict,
        atom_types: list[str],
        low_conf_xyz: str,
    ) -> None:
        min_energy = 1e10
        for ts in trajectory_data:
            if self._save_conformers:
                conformer_file_name = f"conf_{ts}.xyz"
            else:
                # This will get overwrriten each time and deleted at
                # the end.
                conformer_file_name = "temp_conf.xyz"

            self._write_conformer_xyz_file(
                ts=ts,
                ts_data=trajectory_data[ts],
                filename=conformer_file_name,
                atom_types=atom_types,
            )
            mol = mol.with_structure_from_file(conformer_file_name)
            conformer_opt_dir = os.path.join(os.getcwd(), f"conf_{ts}_opt")
            gulp_opt = GulpUFFOptimizer(
                gulp_path=self._gulp_path,
                metal_FF=self._metal_FF,
                output_dir=conformer_opt_dir,
            )
            gulp_opt.assign_FF(mol)
            mol = gulp_opt.optimize(mol=mol)
            energy = gulp_opt.extract_final_energy(
                out_file=os.path.join(conformer_opt_dir, "gulp_opt.ginout"),
            )
            if energy < min_energy:
                min_energy = energy
                # Write out optimised conformer.
                mol.write(low_conf_xyz)

        if not self._save_conformers:
            os.remove("temp_conf.xyz")

    def _save_all_conformers(
        self,
        trajectory_data: dict,
        atom_types: list[str],
    ) -> None:
        for ts in trajectory_data:
            conformer_file_name = f"conf_{ts}.xyz"
            self._write_conformer_xyz_file(
                ts=ts,
                ts_data=trajectory_data[ts],
                filename=conformer_file_name,
                atom_types=atom_types,
            )

    def _calculate_lowest_energy_conformer(
        self,
        trajectory_data: dict[int, dict],
    ) -> int:
        energies = [trajectory_data[ts]["E"] for ts in trajectory_data]
        min_energy = min(energies)
        min_ts = list(trajectory_data.keys())[energies.index(min_energy)]
        return min_ts

    def _save_lowest_energy_conf(
        self,
        mol: stk.Molecule,
        output_xyz: str,
        output_traj: str,
        xyz_traj: str,
        low_conf_xyz: str,
    ) -> None:
        # Convert GULP trajectory file to xyz trajectory.
        (
            atom_types,
            trajectory_data,
            xyz_traj_lines,
        ) = self._convert_traj_to_xyz(
            output_xyz=output_xyz,
            output_traj=output_traj,
        )
        # Write XYZ trajectory file.
        with open(xyz_traj, "w") as f:
            for line in xyz_traj_lines:
                f.write(line)

        # Find lowest energy conformation and output to XYZ.
        if self._opt_conformers:
            self._optimise_all_conformers(
                mol=mol,
                trajectory_data=trajectory_data,
                atom_types=atom_types,
                low_conf_xyz=low_conf_xyz,
            )
        else:
            if self._save_conformers:
                self._save_all_conformers(trajectory_data, atom_types)

            min_ts = self._calculate_lowest_energy_conformer(
                trajectory_data=trajectory_data,
            )

            self._write_conformer_xyz_file(
                ts=min_ts,
                ts_data=trajectory_data[min_ts],
                filename=low_conf_xyz,
                atom_types=atom_types,
            )

    def optimize(self, mol: stk.Molecule) -> stk.Molecule:
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

        in_file = "gulp_MD.gin"
        out_file = "gulp_MD.ginout"
        output_xyz = "gulp_MD_template.xyz"
        output_traj = "gulp_MD.trg"
        xyz_traj = "gulp_MD_traj.xyz"
        low_conf_xyz = "low_energy_conf.xyz"

        metal_atoms = get_metal_atoms(mol)

        try:
            mol.write(output_xyz)
            # Write GULP file.
            self._write_gulp_file(
                mol=mol,
                metal_atoms=metal_atoms,
                in_file=in_file,
                output_traj=output_traj,
            )

            # Run.
            self._run_gulp(in_file, out_file)

            # Get lowest energy conformer from trajectory.
            self._save_lowest_energy_conf(
                mol=mol,
                output_xyz=output_xyz,
                output_traj=output_traj,
                xyz_traj=xyz_traj,
                low_conf_xyz=low_conf_xyz,
            )

            # Update from output.
            mol = mol.with_structure_from_file(low_conf_xyz)

        finally:
            os.chdir(init_dir)

        return mol

    def p_optimize(
        self,
        mol: stk.Molecule,
        unit_cell: UnitCell,
    ) -> tuple[stk.Molecule, UnitCell]:
        """Optimize `mol` and `unit_cell`.

        Parameters
        ----------
            mol:
                The molecule to be optimized.

            unit_cell:
                The unit_cell to be optimized if optimization is periodic.

        Returns:
        -------
            mol:
                The optimized molecule.

            unit_cell:
                The optimized cell.

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

        in_file = "gulp_MD.gin"
        out_file = "gulp_MD.ginout"
        output_xyz = "gulp_MD_template.xyz"
        output_traj = "gulp_MD.trg"
        xyz_traj = "gulp_MD_traj.xyz"
        low_conf_xyz = "low_energy_conf.xyz"
        low_conf_cif = low_conf_xyz.replace("xyz", "cif")

        metal_atoms = get_metal_atoms(mol)

        try:
            mol.write(output_xyz)
            # Write GULP file.
            self._write_gulp_file(
                mol=mol,
                metal_atoms=metal_atoms,
                in_file=in_file,
                output_traj=output_traj,
                unit_cell=unit_cell,
            )

            # Run.
            self._run_gulp(in_file, out_file)

            # Get lowest energy conformer from trajectory.
            self._save_lowest_energy_conf(
                mol=mol,
                output_xyz=output_xyz,
                output_traj=output_traj,
                xyz_traj=xyz_traj,
                low_conf_xyz=low_conf_xyz,
            )

            # Update from output.
            mol = mol.with_structure_from_file(low_conf_xyz)
            unit_cell = unit_cell.with_cell_from_cif(low_conf_cif)

        finally:
            os.chdir(init_dir)

        return mol, unit_cell
