"""
GULP Metal Optimizer
==============

#. :class:`.GulpMetalOptimizer`
#. :class:`.GulpMDMetalOpimizer`

Wrappers for calculators within the :mod:`gulp` code.

Examples
--------

While metal atoms are not required, UFF4MOF is a useful because it
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
        smiles=(
            '[H]C1=NC([H])=C([H])C(C2=C([H])C([H])=C([H])C(C3=C('
            '[H])C([H])=NC([H])=C3[H])=C2[H])=C1[H]'
        ),
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

    gulp_opt = stko.GulpMetalOptimizer(
        gulp_path='path/to/gulp',
        metal_FF={46: 'Pd4+2'},
        conjugate_gradient=True
    )

    # Assign the force field.
    gulp_opt.assign_FF(cage)
    # Run optimization.
    cage = gulp_opt.optimize(mol=cage)

Conformer searching is often useful, so we have provided an interface
to MD simulations using GULP and UFF4MOF. A conformer search can be
run at high temperature, where N conformers are extracted at constant
intervals throughtout the simulation and optimized using UFF4MOF. The
lowest energy conformer is returned. After these MD steps, it is
crucial to reoptimize the resultant structure using a better forcefield
or a more robust method!

.. code-block:: python

    gulp_MD = stko.GulpMDMetalOptimizer(
        gulp_path='path/to/gulp',
        metal_FF={46: 'Pd4+2'},
        temperature=700,
        N_conformers=10,
        opt_conformers=True,
    )
    gulp_MD.assign_FF(cage)
    cage = gulp_MD.optimize(cage)

"""

import logging
import re
import uuid
import os
import shutil
import subprocess as sp
from rdkit.Chem import AllChem as rdkit

from .optimizers import Optimizer
from ..utilities import (
    has_h_atom,
    has_metal_atom,
    get_metal_atoms,
    get_metal_bonds,
    to_rdkit_mol_without_metals
)


logger = logging.getLogger(__name__)


class UFFTyperError(Exception):
    ...


class GulpMetalOptimizer(Optimizer):
    """
    Applies forcefield optimizers that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run an optimisation using the
    UFF4MOF. This forcefield requires some explicit metal atom
    definitions, which are determined by the user.

    """

    def __init__(
        self,
        gulp_path,
        metal_FF,
        conjugate_gradient=False,
        output_dir=None,
    ):
        """
        Initialize a :class:`GulpMetalOptimizer` instance.

        Parameters
        ----------
        gulp_path : :class:`str`
            Path to GULP executable.

        metal_FF : :class:`dict`
            Dictionary with metal atom forcefield assignments.
            Key: :class:`int` : atomic number.
            Value: :class:`str` : UFF4MOF forcefield type.

        conjugate_gradient : :class:``, optional
            ``True`` to use Conjugate Graditent method.
            Defaults to ``False``

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        """
        self._gulp_path = gulp_path
        self._metal_FF = metal_FF
        self._conjugate_gradient = conjugate_gradient
        self._output_dir = output_dir

    def _add_atom_charge_flags(self, atom, atomkey):
        """
        Add atom charge flags for forcefield.

        Code inspired by:
        https://github.com/rdkit/rdkit/blob/master/Code/
        GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        total_valence = rdkit.Atom.GetTotalValence(atom)
        atnum = int(atom.GetAtomicNum())

        # Go through element cases.
        # Mg.
        if atnum == 12:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Al.
        elif atnum == 13:
            if total_valence != 3:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Si.
        elif atnum == 14:
            if total_valence != 4:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # P.
        elif atnum == 15:
            if total_valence == 3:
                atomkey += '+3'
            elif total_valence == 5:
                atomkey += '+5'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # S.
        elif atnum == 16:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if hybrid != rdkit.HybridizationType.SP2:
                if total_valence == 2:
                    atomkey += '+2'
                elif total_valence == 4:
                    atomkey += '+4'
                elif total_valence == 6:
                    atomkey += '+6'
                else:
                    raise UFFTyperError(
                        f"UFFTYPER: Unrecognized charge state for "
                        f"atom: {atom.GetIdx}"
                    )
        # Zn.
        elif atnum == 30:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Ga.
        elif atnum == 31:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # As.
        elif atnum == 33:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Se.
        elif atnum == 34:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Cd.
        elif atnum == 48:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # In.
        elif atnum == 49:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )

        # Sb.
        elif atnum == 51:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Te.
        elif atnum == 52:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Hg.
        elif atnum == 80:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Tl.
        elif atnum == 81:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Pb.
        elif atnum == 82:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Bi.
        elif atnum == 83:
            if total_valence == 3:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Po.
        elif atnum == 84:
            if total_valence == 2:
                atomkey += '+2'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        # Lanthanides.
        elif atnum >= 57 and atnum <= 71:
            if total_valence == 6:
                atomkey += '+3'
            else:
                raise UFFTyperError(
                    f"UFFTYPER: Unrecognized charge state for atom: "
                    f"{atom.GetIdx}"
                )
        return atomkey

    def _get_atom_label(self, atom):
        """
        Get FF atom label.

        Code inspired by:
        https://github.com/rdkit/rdkit/blob/master/Code/
        GraphMol/ForceFieldHelpers/UFF/AtomTyper.cpp

        """
        atnum = int(atom.GetAtomicNum())
        atomkey = atom.GetSymbol()
        if len(atomkey) == 1:
            atomkey += '_'

        table = rdkit.GetPeriodicTable()

        chk1 = (
            rdkit.PeriodicTable.GetDefaultValence(table, atnum) == -1
        )
        chk2 = (rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 1)
        chk3 = (rdkit.PeriodicTable.GetNOuterElecs(table, atnum) != 7)
        chk4 = chk2 and chk3
        if chk1 or chk4:
            hybrid = rdkit.Atom.GetHybridization(atom)
            if atnum == 84:
                atomkey += '3'
                if hybrid != rdkit.HybridizationType.SP3:
                    raise UFFTyperError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            elif atnum == 80:
                atomkey += '1'
                if hybrid != rdkit.HybridizationType.SP:
                    raise UFFTyperError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
            else:
                if hybrid == rdkit.HybridizationType.SP:
                    atomkey += '1'
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
                    chk4a = (chk1a or chk2a)
                    if chk4a and chk3a:
                        atomkey += 'R'
                    else:
                        atomkey += '2'
                elif hybrid == rdkit.HybridizationType.SP3:
                    atomkey += '3'
                elif hybrid == rdkit.HybridizationType.SP3D:
                    atomkey += '5'
                elif hybrid == rdkit.HybridizationType.SP3D2:
                    atomkey += '6'
                else:
                    raise UFFTyperError(
                        f"UFFTYPER: Unrecognized hybridization for"
                        f" atom: {atom.GetIdx}"
                    )
        atomkey = self._add_atom_charge_flags(atom, atomkey)
        return atomkey

    def _type_translator(self):
        type_translator = {}
        types = set([self.atom_labels[i][0] for i in self.atom_labels])
        for t in types:
            if not t[1].isalpha():
                symb = t[0]
            else:
                symb = t[0:2]
            for i in range(1, 100):
                name = f'{symb}{i}'
                if name in type_translator.values():
                    continue
                else:
                    type_translator[t] = name
                    break

        return type_translator

    def _position_section(self, mol, type_translator):
        position_section = '\ncartesian\n'
        for atom in mol.get_atoms():
            atom_type = type_translator[self.atom_labels[
                atom.get_id()
            ][0]]
            position = mol.get_centroid(atom_ids=atom.get_id())
            posi_string = (
                f'{atom_type} core {round(position[0], 5)} '
                f'{round(position[1], 5)} {round(position[2], 5)}\n'
            )
            position_section += posi_string

        return position_section

    def _bond_section(self, mol, metal_atoms):
        bond_section = '\n'
        for bond in mol.get_bonds():
            atom_types = [
                self.atom_labels[i.get_id()][0]
                for i in [bond.get_atom1(), bond.get_atom2()]
            ]

            # Set bond orders.
            if has_h_atom(bond):
                # H has bond order of 1.
                bond_type = ''
            elif has_metal_atom(bond, metal_atoms):
                bond_type = 'half'
            elif '_R' in atom_types[0] and '_R' in atom_types[1]:
                bond_type = 'resonant'
            elif bond.get_order() == 1:
                bond_type = ''
            elif bond.get_order() == 2:
                bond_type = 'double'
            elif bond.get_order() == 3:
                bond_type = 'triple'

            string = (
                f'connect {bond.get_atom1().get_id()+1} '
                f'{bond.get_atom2().get_id()+1} {bond_type}'
            )
            bond_section += string+'\n'

        return bond_section

    def _species_section(self, type_translator):
        species_section = '\nspecies\n'
        for spec in type_translator:
            name = type_translator[spec]
            species_section += f'{name} {spec}\n'

        return species_section

    def _write_gulp_file(self, mol, metal_atoms, in_file, output_xyz):

        type_translator = self._type_translator()

        if self._conjugate_gradient:
            top_line = (
                'opti conj unit conv noautobond fix molmec cartesian\n'
            )
        else:
            top_line = 'opti conv noautobond fix molmec cartesian\n'

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = '\nlibrary uff4mof.lib\n'

        output_section = (
            '\n'
            'terse inout potentials\n'
            'terse inout cell\n'
            'terse in structure\n'
            'terse inout derivatives\n'
            f'output xyz {output_xyz}\n'
            # 'output movie xyz steps_.xyz\n'
        )

        with open(in_file, 'w') as f:
            f.write(top_line)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(output_section)

    def assign_FF(self, mol):
        """
        Assign forcefield types to molecule.

        Parameters
        ----------

        Returns
        -------

        None : :class:`NoneType`

        """
        metal_atoms = get_metal_atoms(mol)
        metal_ids = [i.get_id() for i in metal_atoms]
        metal_bonds, _ = get_metal_bonds(mol, metal_atoms)
        edit_mol = to_rdkit_mol_without_metals(
            mol=mol,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )

        # Get forcefield parameters.
        rdkit.SanitizeMol(edit_mol)
        self.atom_labels = {}

        for i in range(edit_mol.GetNumAtoms()):
            if i in metal_ids:
                self.atom_labels[i] = [None, 'metal', None]
            else:
                atom = edit_mol.GetAtomWithIdx(i)
                atom_label = self._get_atom_label(atom)
                self.atom_labels[i] = [atom_label, None, None]

        # Write UFF4MOF specific forcefield parameters.
        # Metals.
        for atomid in self.atom_labels:
            if self.atom_labels[atomid][1] == 'metal':
                atom, = mol.get_atoms(atomid)
                atom_no = atom.get_atomic_number()
                self.atom_labels[atomid][0] = self._metal_FF[atom_no]

    def _run_gulp(self, in_file, out_file):
        cmd = f'{self._gulp_path} < {in_file}'
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

    def _move_generated_files(self, files):
        if not os.path.exists(self._output_dir):
            os.mkdir(self._output_dir)

        for file in files:
            if os.path.exists(file):
                os.rename(file, f'{self._output_dir}/{file}')

    def extract_final_energy(self, file):
        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        with open(file, 'r') as f:
            for line in f.readlines():
                if 'Final energy =' in line:
                    string = nums.search(line.rstrip()).group(0)
                    return float(string)

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

        in_file = 'gulp_opt.gin'
        out_file = 'gulp_opt.ginout'
        output_xyz = 'gulp_opt.xyz'

        metal_atoms = get_metal_atoms(mol)

        # Write GULP file.
        self._write_gulp_file(
            mol=mol,
            metal_atoms=metal_atoms,
            in_file=in_file,
            output_xyz=output_xyz
        )

        # Run.
        self._run_gulp(in_file, out_file)

        # Update from output.
        mol = mol.with_structure_from_file(output_xyz)

        # Move files.
        self._move_generated_files(
            files=[in_file, out_file, output_xyz]
        )

        return mol


class GulpMDMetalOptimizer(GulpMetalOptimizer):
    """
    Applies forcefield MD that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run a MD run using the UFF4MOF.
    This forcefield requires some explicit metal atom definitions,
    which are determined by the user.

    """

    def __init__(
        self,
        gulp_path,
        metal_FF,
        output_dir=None,
        integrator='stochastic',
        ensemble='nvt',
        temperature=300,
        equilbration=1.0,
        production=10.0,
        timestep=1.0,
        N_conformers=10,
        opt_conformers=True,
        save_conformers=False,
    ):
        """
        Initialize a :class:`GulpMetalOptimizer` instance.

        Parameters
        ----------
        gulp_path : :class:`str`
            Path to GULP executable.

        metal_FF : :class:`dict`
            Dictionary with metal atom forcefield assignments.
            Key: :class:`int` : atomic number.
            Value: :class:`str` : UFF4MOF forcefield type.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        integrator : :class:`str`, optional
            Integrator for GULP to use.
            Defaults to 'stochastic'.

        ensemble : :class:`str`, optional
            Ensemble for GULP to use.
            Defaults to 'nvt'.

        temperature : :class:`float`, optional
            Temperature to run simulation at in Kelvin.
            Defaults to 300.

        equilbration : :class:`float`, optional
            Time spent equilibrating system in ps.
            Defaults to 1.0.

        production : :class:`float`, optional
            Time spent running production simulation of system in ps.
            Defaults to 10.0.

        timestep : :class:`float`, optional
            Timestep of simulation in fs.
            Defaults to 1.0.

        N_conformers : :class:`int`, optional
            Number of conformers to sample.
            Defaults to 10.

        opt_conformers : :class:`bool`, optional
            Whether or not to optimise each conformer using UFF4MOF.
            Defaults to ``True``.

        save_conformers : :class:`bool`, optional
            Whether or not to save to file each conformer.
            Defaults to ``False``.

        """
        self._gulp_path = gulp_path
        self._metal_FF = metal_FF
        self._output_dir = output_dir
        self._integrator = integrator
        self._ensemble = ensemble
        self._temperature = temperature
        self._equilbration = float(equilbration)
        self._production = float(production)
        self._timestep = timestep
        self._N_conformers = N_conformers
        samples = float(self._production) / float(self._N_conformers)
        self._sample = samples
        self._write = samples
        self._opt_conformers = opt_conformers
        self._save_conformers = save_conformers

    def _write_gulp_file(
        self,
        mol,
        metal_atoms,
        in_file,
        output_traj
    ):

        type_translator = self._type_translator()

        top_line = 'md conv noautobond fix molmec cartesian\n'

        position_section = self._position_section(mol, type_translator)
        bond_section = self._bond_section(mol, metal_atoms)
        species_section = self._species_section(type_translator)

        library = '\nlibrary uff4mof.lib\n'

        md_section = (
            '\n'
            "mdmaxtemp 100000000\n"
            f"integrator {self._integrator}\n"
            f"ensemble {self._ensemble}\n"
            f"temperature {self._temperature}\n"
            f"equilbration {self._equilbration} ps\n"
            f"production {self._production} ps\n"
            f"timestep {self._timestep} fs\n"
            f"sample {self._sample} ps\n"
            f"write {self._write} ps\n"
            '\n'
        )

        output_section = (
            '\n'
            'terse inout potentials\n'
            'terse inout cell\n'
            'terse in structure\n'
            'terse inout derivatives\n'
            f'output trajectory ascii {output_traj}\n'
        )

        with open(in_file, 'w') as f:
            f.write(top_line)
            f.write(position_section)
            f.write(bond_section)
            f.write(species_section)
            f.write(library)
            f.write(md_section)
            f.write(output_section)

    def _convert_traj_to_xyz(self, output_xyz, output_traj, xyz_traj):

        # Get atom types from an existing xyz file.
        atom_types = []
        with open(output_xyz, 'r') as f:
            for line in f.readlines()[2:]:
                atom_types.append(line.rstrip().split(' ')[0])

        # Read in lines from trajectory file.
        with open(output_traj, 'r') as f:
            lines = f.readlines()

        # Split file using strings.
        id = 0
        s_times = {}
        s_coords = {}
        s_vels = {}
        s_derivs = {}
        s_sites = {}
        coords = []
        vels = []
        derivs = []
        sites = []
        switch = None
        for line in lines:
            li = line.rstrip()
            if '#  Time/KE/E/T' in li:
                switch = 'times'
                s_coords[id] = coords
                coords = []
                s_vels[id] = vels
                vels = []
                s_derivs[id] = derivs
                derivs = []
                s_sites[id] = sites
                sites = []
                id += 1
            elif '#  Coordinates' in li:
                switch = 'coords'
            elif '#  Velocities' in li:
                switch = 'vels'
            elif '#  Derivatives' in li:
                switch = 'derivs'
            elif '#  Site energies' in li:
                switch = 'sites'
            elif switch == 'coords':
                coords.append([i for i in li.split(' ') if i])
            elif switch == 'vels':
                vels.append([i for i in li.split(' ') if i])
            elif switch == 'derivs':
                derivs.append([i for i in li.split(' ') if i])
            elif switch == 'sites':
                sites.append([i for i in li.split(' ') if i])
            elif switch == 'times':
                s_times[id] = [i for i in li.split(' ') if i]
            elif switch is None:
                pass
        # Add final timestep.
        s_coords[id] = coords
        coords = []
        s_vels[id] = vels
        vels = []
        s_derivs[id] = derivs
        derivs = []
        s_sites[id] = sites
        sites = []

        ids = []
        tt = []
        pot_energies = []
        new_lines = []
        for id in s_times:
            times = s_times[id]
            ids.append(id)
            tt.append(float(times[0]))
            pot_energies.append(float(times[2]))
            coords = s_coords[id]
            sites = s_sites[id]
            xyz_string = (
                f'{len(coords)}\n'
                f'{times[0]},{times[1]},{times[2]},{times[3]}\n'
            )

            for i, coord in enumerate(coords):
                site = sites[i][0]
                xyz_string += (
                    f'{atom_types[i]} {coord[0]} {coord[1]} '
                    f'{coord[2]} {site}\n'
                )

            new_lines.append(xyz_string)

        with open(xyz_traj, 'w') as f:
            for line in new_lines:
                f.write(line)

        return atom_types, ids, tt, pot_energies, s_times, s_coords

    def _write_conformer_xyz_file(
        self,
        id,
        filename,
        s_times,
        s_coords,
        atom_types
    ):
        times = s_times[id]
        coords = s_coords[id]
        xyz_string = (
            f'{len(coords)}\n'
            f'{times[0]},{times[1]},{times[2]},{times[3]}\n'
        )
        for i, coord in enumerate(coords):
            xyz_string += (
                f'{atom_types[i]} {coord[0]} {coord[1]} {coord[2]}\n'
            )
        with open(filename, 'w') as f:
            f.write(xyz_string)

    def _save_lowest_energy_conf(
        self,
        mol,
        output_xyz,
        output_traj,
        xyz_traj,
        low_conf_xyz
    ):

        # Convert GULP trajectory file to xyz trajectory.
        results = self._convert_traj_to_xyz(
            output_xyz,
            output_traj,
            xyz_traj
        )
        atom_types, ids, tt, pot_energies, s_times, s_coords = results

        # Find lowest energy conformation and output to XYZ.
        if self._opt_conformers:
            # Optimise all conformers.
            min_energy = 1E10
            for id in ids:
                if self._save_conformers:
                    conformer_file_name = os.path.join(
                        self._output_dir,
                        f'conf_{id}.xyz'
                    )
                else:
                    # This will get overwrriten each time.
                    conformer_file_name = os.path.join(
                        self._output_dir,
                        'temp_conf.xyz'
                    )
                self._write_conformer_xyz_file(
                    id=id,
                    filename=conformer_file_name,
                    s_times=s_times,
                    s_coords=s_coords,
                    atom_types=atom_types
                )
                mol = mol.with_structure_from_file(conformer_file_name)
                gulp_opt = GulpMetalOptimizer(
                    gulp_path=self._gulp_path,
                    metal_FF=self._metal_FF,
                    output_dir=self._output_dir
                )
                gulp_opt.assign_FF(mol)
                mol = gulp_opt.optimize(mol=mol)
                energy = gulp_opt.extract_final_energy(
                    file=os.path.join(
                        self._output_dir,
                        'gulp_opt.ginout'
                    )
                )

                if energy < min_energy:
                    min_energy = energy
                    self._write_conformer_xyz_file(
                        id=id,
                        filename=low_conf_xyz,
                        s_times=s_times,
                        s_coords=s_coords,
                        atom_types=atom_types
                    )
        else:
            if self._save_conformers:
                for id in ids:
                    conformer_file_name = os.path.join(
                        self._output_dir,
                        f'conf_{id}.xyz'
                    )
                    self._write_conformer_xyz_file(
                        id=id,
                        filename=conformer_file_name,
                        s_times=s_times,
                        s_coords=s_coords,
                        atom_types=atom_types
                    )

            min_energy = min(pot_energies)
            min_id = ids[pot_energies.index(min_energy)]
            self._write_conformer_xyz_file(
                id=min_id,
                filename=low_conf_xyz,
                s_times=s_times,
                s_coords=s_coords,
                atom_types=atom_types
            )

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

        in_file = 'gulp_MD.gin'
        out_file = 'gulp_MD.ginout'
        output_xyz = 'gulp_MD_template.xyz'
        mol.write(output_xyz)
        output_traj = 'gulp_MD.trg'
        xyz_traj = 'gulp_MD_traj.xyz'
        low_conf_xyz = 'low_energy_conf.xyz'

        metal_atoms = get_metal_atoms(mol)

        # Write GULP file.
        self._write_gulp_file(
            mol=mol,
            metal_atoms=metal_atoms,
            in_file=in_file,
            output_traj=output_traj
        )

        # Run.
        self._run_gulp(in_file, out_file)

        # Get lowest energy conformer from trajectory.
        self._save_lowest_energy_conf(
            mol=mol,
            output_xyz=output_xyz,
            output_traj=output_traj,
            xyz_traj=xyz_traj,
            low_conf_xyz=low_conf_xyz
        )

        # Update from output.
        mol = mol.with_structure_from_file(low_conf_xyz)

        # Move files.
        self._move_generated_files(
            files=[
                in_file, out_file, output_xyz,
                output_traj, xyz_traj, low_conf_xyz,
                'temp_conf.xyz'
            ]
        )

        return mol
