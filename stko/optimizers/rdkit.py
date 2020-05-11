"""
RDKit Optimizers
================

#. :class:`.MMFF`
#. :class:`.UFF`
#. :class:`.ETKDG`
#. :class:`.MetalOptimizer`

Wrappers for optimizers within the :mod:`rdkit` code.

.. code-block:: python

    import stk
    import stko

    # Optimizers work on stk.Molecule.
    polymer = stk.ConstructedMolecule(
        building_blocks=[mol],
        topology_graph=stk.polymer.Linear('A', [0], n=3)
    )
    etkdg = stk.ETKDG()
    mol = etkdg.optimize(polymer)

"""

import logging
import rdkit.Chem.AllChem as rdkit

from .optimizers import Optimizer


logger = logging.getLogger(__name__)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        mmff = stk.MMFF()
        mol = mmff.optimize(mol)

    """

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(rdkit_mol)
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        uff = stk.UFF()
        mol = uff.optimize(mol)

    """

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

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(rdkit_mol)
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )

        return mol


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        etkdg = stk.ETKDG()
        mol = etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12):
        """
        Initialize a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        """

        self._random_seed = random_seed

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

        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        mol = mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )

        return mol


class MetalOptimizer(Optimizer):
    """
    Applies forcefield optimizers that can handle metal centres.

    Notes
    -----
    By default, :meth:`optimize` will run a restricted optimization
    using constraints and the UFF. To implement this, metal atoms are
    replaced by noninteracting H atoms, and constraints are applied
    to maintain the metal centre geometry.

    Restrictions are applied to the ligand with respect to its input
    structure. So if that is poorly optimised, then the output will be
    also.

    Examples
    --------

    :class:`MetalOptimizer` allows for the restricted optimization of
    :class:`ConstructedMolecule` containing metals. Note that this
    optimizer algorithm is not very robust to large bonds and may fail.

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
        cage1 = stk.ConstructedMolecule(
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

        # Define an optimizer, with ligand bonder atoms placed 1.9
        # Angstrom from the metal centres. A weak force_constant is
        # used here to not overshadow all other forces. A slow
        # optimisation of the ligand - bonder bonds is done over
        # res_steps, where each step is a very short RDKit UFF
        # optimisation with constraints.
        # All ligand bonds and angles have restricitons applied based
        # on the input structure.
        optimizer = stk.MetalOptimizer(
            scale=1.9,
            force_constant=1.0e2,
            rel_distance=0.9,
            res_steps=50,
            restrict_all_bonds=True,
            restrict_orientation=True
        )

        # Optimize.
        cage1 = optimizer.optimize(mol=cage1)

    """

    def __init__(
        self,
        metal_binder_distance,
        metal_binder_fc,
        binder_ligand_fc,
        ignore_vdw,
        rel_distance,
        res_steps,
        max_iterations,
        do_long_opt,
        restrict_bonds=False,
        restrict_angles=False,
        restrict_orientation=False
    ):
        """
        Initialize a :class:`MetalOptimizer` instance.

        Parameters
        ----------

        force_constant : :class:`float`
            Force constant to use for restricted metal-ligand bonds.
            All other force constants are hard-coded/standard values.

        rel_distance : :class:`float`
            Set the relative distance to optimise the metal-ligand
            bonds to.

        res_steps : :class:`int`
            Number of optimisation steps to run. Each optimisation step
            is very short to ensure the structure does not over-react
            to the extreme forces on it.

        restrict_all_bonds : :class:`bool`, optional
            `True` to restrict all bonds except for ligand-FG bonds.

        restrict_orientation : :class:`bool`, optional
            `True` to restrict metal complex FG angles relative to
            topology centre of mass.

        """
        self._metal_binder_distance = metal_binder_distance
        self._metal_binder_fc = metal_binder_fc
        self._binder_ligand_fc = binder_ligand_fc
        self._ignore_vdw = ignore_vdw
        self._rel_distance = rel_distance
        self._restrict_bonds = restrict_bonds
        self._restrict_angles = restrict_angles
        self._restrict_orientation = restrict_orientation
        self._res_steps = res_steps
        self._max_iterations = max_iterations
        self._do_long_opt = do_long_opt

    def get_input_constraints(
        self,
        mol,
        ids_to_metals,
        metal_atoms,
        include_bonders=False
    ):
        """
        Get a series of constraint definitions based on mol.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        ids_to_metals : :class:`.list` of :class:`int`
            List of atom ids of atoms bonded to metals.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        include_bonders : :class:`bool`, optional
            Whether to constrain the angles of atoms bonded to the
            metal atoms. Defaults to `False`

        Returns
        -------
        constraints : :class:`dict`
            Dictionary of all constraints of the form:
                bonds: {stk.bond: {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'type': 'bond',
                    'fc': :class:`float`
                }}
                angles: {(stk.bond2, stk.bond): {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'idx3': :class:`int`,
                    'type': 'angle',
                    'angle': :class:`float`,
                    'fc': :class:`float`
                }}
                torsions: {(stk.bond, stk.bond, stk.bond): {
                    'idx1': :class:`int`,
                    'idx2': :class:`int`,
                    'idx3': :class:`int`,
                    'idx4': :class:`int`,
                    'type': 'torsion',
                    'torsion': :class:`float`,
                    'fc': :class:`float`
                }}

        """
        constraints = {}
        # Bond constraints. Set to have a force constant of 1E2.
        for bond in mol.bonds:
            idx1 = bond.atom1.id
            idx2 = bond.atom2.id
            if idx1 in ids_to_metals or idx2 in ids_to_metals:
                continue
            # Do not restrict H atom distances.
            if self.has_h(bond):
                continue
            if self.has_M(bond, metal_atoms):
                continue
            length = mol.get_atom_distance(idx1, idx2)
            constraints[bond] = {
                'idx1': idx1,
                'idx2': idx2,
                'type': 'bond',
                'length': length,
                'fc': 1.0e2
            }

        # Angle constraints
        # Add angle constraints to angles containing H atoms.
        # Set to have a force constant of 1E1, unless the angle
        # includes H atoms, then it has a force constant of 1E2.
        for bonds in combinations(mol.bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) > 3:
                continue
            # Only want angles with at least 1 X-H bond.
            # if not self.has_h(bond1) and not self.has_h(bond2):
            #     continue
            # Do not want any metal containing bonds.
            if any([self.has_M(i, metal_atoms) for i in bonds]):
                continue
            for atom in pres_atoms:
                if atom in bond1_atoms and atom in bond2_atoms:
                    idx2 = atom.id
                elif atom in bond1_atoms:
                    idx1 = atom.id
                elif atom in bond2_atoms:
                    idx3 = atom.id
            if not include_bonders:
                # Do not restrict angles of bonds to atoms bonded to
                # metals.
                bonded_to_metals = [
                    idx1 in ids_to_metals,
                    idx2 in ids_to_metals,
                    idx3 in ids_to_metals
                ]
                if any(bonded_to_metals):
                    continue
            pos1 = [
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            if self.has_h(bond1) or self.has_h(bond2):
                # Increased force constant for angles with H atoms.
                constraints[(bond1, bond2)] = {
                    'idx1': idx1,
                    'idx2': idx2,
                    'idx3': idx3,
                    'type': 'angle',
                    'angle': np.degrees(angle),
                    'fc': 1.0e2
                }
            else:
                constraints[(bond1, bond2)] = {
                    'idx1': idx1,
                    'idx2': idx2,
                    'idx3': idx3,
                    'type': 'angle',
                    'angle': np.degrees(angle),
                    'fc': 1.0e1
                }

        return constraints

    def _restricted_optimization(
        self,
        ff,
        mol,
        edit_mol,
        metal_atoms,
        metal_bonds,
        ids_to_metals,
        input_constraints=None
    ):
        """
        Optimize `mol` with restrictions on metal-ligand bonds.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        ids_to_metals : :class:`.list` of :class:`int`
            List of atom ids of atoms bonded to metals.

        input_constraints : :class:`dict`
            Dictionary of constraints to apply.

        Returns
        -------
        None : :class:`NoneType`

        """

        # For bonds between ligand bonders and the rest of the ligand,
        # a weak force constant is applied to minimize to rel_distance.
        # This is the slow relaxation of the high-force bonds.
        if self._rel_distance is not None:
            for bond in mol.bonds:
                idx1 = bond.atom1.id
                idx2 = bond.atom2.id
                if idx1 in ids_to_metals or idx2 in ids_to_metals:
                    # Do not restrict H atom distances.
                    if self.has_h(bond):
                        continue
                    if self.has_M(bond, metal_atoms):
                        continue
                    # Add distance constraints in place of metal bonds.
                    distance = mol.get_atom_distance(idx1, idx2)
                    ff.UFFAddDistanceConstraint(
                        idx1=idx1,
                        idx2=idx2,
                        relative=False,
                        minLen=self._rel_distance*distance,
                        maxLen=self._rel_distance*distance,
                        forceConstant=self._binder_ligand_fc
                    )

        # Perform UFF optimization with rdkit.
        ff.Minimize(maxIts=self._max_iterations)

        # Update stk molecule from optimized molecule. This should
        # only modify atom positions, which means metal atoms will be
        # reinstated.
        mol.update_from_rdkit_mol(edit_mol)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Find all metal atoms and atoms they are bonded to.
        metal_atoms = self.get_metal_atoms(mol)
        metal_bonds, ids_to_metals = self.get_metal_bonds(
            mol,
            metal_atoms
        )

        input_constraints = self.get_input_constraints(
            mol,
            ids_to_metals,
            metal_atoms,
            include_bonders=False
        )

        # Perform a forcefield optimisation that
        # only optimises non metal atoms that are not bonded to the
        # metal.

        # Write rdkit molecule with metal atoms and bonds deleted.
        edit_mol = self.to_rdkit_mol_no_metals(
            mol,
            metal_atoms=metal_atoms,
            metal_bonds=metal_bonds
        )

        # Non-bonded interaction interactions need to be turned on
        # because at the point of intialisation, the molecule are
        # technically separate (not bonded) and are treated as
        # fragments.
        rdkit.SanitizeMol(edit_mol)
        ff = rdkit.UFFGetMoleculeForceField(
            edit_mol,
            ignoreInterfragInteractions=self._ignore_vdw
        )

        # Constrain selected bonds, angles and dihedrals and metal
        # atoms.
        # Constrain the metal centre.
        self.apply_metal_centre_constraints(
            mol,
            ff,
            metal_bonds,
            metal_atoms
        )

        # Add angular constraints that enforce relative orientation
        # between metal complexes + the topology centre of mass.
        if self._restrict_orientation:
            self.apply_orientation_restriction(
                ff,
                mol,
                metal_bonds,
                metal_atoms
            )

        # Constrain all bonds and angles based on input structure
        # except for:
        # (1) bonds including metals
        # (2) bonds including atoms bonded to metals
        for constraint in input_constraints:
            const = input_constraints[constraint]
            if const['type'] == 'bond' and self._restrict_bonds:
                # Add distance constraints in place of metal bonds.
                ff.UFFAddDistanceConstraint(
                    idx1=const['idx1'],
                    idx2=const['idx2'],
                    relative=False,
                    minLen=const['length'],
                    maxLen=const['length'],
                    forceConstant=const['fc']
                )
            elif const['type'] == 'angle' and self._restrict_angles:
                ff.UFFAddAngleConstraint(
                    idx1=const['idx1'],
                    idx2=const['idx2'],
                    idx3=const['idx3'],
                    relative=False,
                    minAngleDeg=const['angle'],
                    maxAngleDeg=const['angle'],
                    forceConstant=const['fc']
                )

        # Optimisation with UFF in RDKit. This method uses constraints
        # on the metal centre to attempt to enforce the metal geometry
        # described by the metal topology.
        # For RDKit, the optimisation is done in loops, where the
        # metal-ligand and adjacent bonds are slowly relaxed to normal
        # values.
        # The rest of the ligands are constrained to the input value.
        for i in range(self._res_steps):
            self._restricted_optimization(
                mol=mol,
                edit_mol=edit_mol,
                ff=ff,
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                input_constraints=input_constraints
            )

        # Finish with one long optimisation.
        if self._do_long_opt:
            self._max_iterations = 500
            self._restricted_optimization(
                mol=mol,
                edit_mol=edit_mol,
                ff=ff,
                metal_atoms=metal_atoms,
                metal_bonds=metal_bonds,
                ids_to_metals=ids_to_metals,
                input_constraints=input_constraints
            )

    def apply_metal_centre_constraints(
        self,
        mol,
        ff,
        metal_bonds,
        metal_atoms
    ):
        """
        Applies UFF metal centre constraints.

        Parameters
        ----------
        ff : :class:`rdkit.ForceField`
            Forcefield to apply constraints to. Generally use UFF.

        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Add a very weak force constraint on all metal-metal
        # distances.
        #
        # for atoms in combinations(metal_atoms, r=2):
        #     ff.UFFAddDistanceConstraint(
        #         idx1=atoms[0].id,
        #         idx2=atoms[1].id,
        #         relative=True,
        #         minLen=0.8,
        #         maxLen=0.8,
        #         forceConstant=0.25e1
        #     )

        # Add constraints to UFF to hold metal geometry in place.
        for bond in metal_bonds:
            idx1 = bond.atom1.id
            idx2 = bond.atom2.id
            # Add distance constraints in place of metal bonds.
            # Target distance set to a given metal_binder_distance.
            ff.UFFAddDistanceConstraint(
                idx1=idx1,
                idx2=idx2,
                relative=False,
                minLen=self._metal_binder_distance,
                maxLen=self._metal_binder_distance,
                forceConstant=self._metal_binder_fc
            )

        # Also implement angular constraints to all atoms in the
        # metal complex.
        for bonds in combinations(metal_bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) > 3:
                continue
            for atom in pres_atoms:
                if atom in bond1_atoms and atom in bond2_atoms:
                    idx2 = atom.id
                elif atom in bond1_atoms:
                    idx1 = atom.id
                elif atom in bond2_atoms:
                    idx3 = atom.id
            pos1 = [
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle),
                maxAngleDeg=np.degrees(angle),
                forceConstant=1.0e5
            )

    def apply_orientation_restriction(
        self,
        ff,
        mol,
        metal_bonds,
        metal_atoms
    ):
        """
        Applies UFF relative orientation restrcitions.

        Parameters
        ----------
        ff : :class:`rdkit.ForceField`
            Forcefield to apply constraints to. Generally use UFF.

        mol : :class:`.Molecule`
            The molecule to be optimized.

        metal_bonds : :class:`.list` of :class:`stk.Bond`
            List of bonds including metal atoms.

        metal_atoms : :class:`.list` of :class:`stk.Atom`
            List of metal atoms.

        Returns
        -------
        None : :class:`NoneType`

        """
        # Add a fixed point.
        COM = mol.get_center_of_mass()
        nidx = ff.AddExtraPoint(
            x=COM[0],
            y=COM[1],
            z=COM[2],
            fixed=True
        )
        ff.Initialize()
        for bond in metal_bonds:
            if bond.atom1 in metal_atoms:
                idx2 = bond.atom1.id
                idx1 = bond.atom2.id
            elif bond.atom2 in metal_atoms:
                idx1 = bond.atom1.id
                idx2 = bond.atom2.id
            pos1 = [
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = COM
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=nidx-1,
                relative=False,
                minAngleDeg=np.degrees(angle)-2,
                maxAngleDeg=np.degrees(angle)+2,
                forceConstant=1.0e4
            )

        # Apply an angular constraint on the binder-metal-metal atoms
        # to maintain the metal centres relative orientation.
        for bonds in combinations(metal_bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.atom1, bond1.atom2]
            bond2_atoms = [bond2.atom1, bond2.atom2]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) < 4:
                continue

            if bond1.atom1 in metal_atoms:
                idx1 = bond1.atom1.id
                idx2 = bond1.atom2.id
            else:
                idx1 = bond1.atom2.id
                idx2 = bond1.atom1.id

            if bond2.atom1 in metal_atoms:
                idx3 = bond2.atom1.id
            else:
                idx3 = bond2.atom2.id

            pos1 = [
                i for i in mol.get_atom_positions(atom_ids=[idx1])
            ][0]
            pos2 = [
                i for i in mol.get_atom_positions(atom_ids=[idx2])
            ][0]
            pos3 = [
                i for i in mol.get_atom_positions(atom_ids=[idx3])
            ][0]
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle),
                maxAngleDeg=np.degrees(angle),
                forceConstant=1.0e4
            )
