"""
Collapser Optimizer
===================

#. :class:`.Collapser`

Optimizer for collapsing enlarged topologies.

"""

import logging
from itertools import combinations
import numpy as np
import uuid
import os
import shutil

from .optimizers import Optimizer
from ..utilities import get_atom_distance


logger = logging.getLogger(__name__)


class Collapser(Optimizer):
    """
    Collapse stk.ConstructedMolecule to decrease enlarged bonds.

    This optimizer aims to bring extended bonds closer together for
    further optimisation.

    .. code-block:: python

        import stk
        import stko

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage1 = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
        )

        # Perform collapser optimisation.
        output_dir = f'cage_opt_{cage_name}_coll'
        optimizer = stko.Collapser(
            output_dir=output_dir,
            step_size=0.05,
            distance_cut=2.0,
            scale_steps=True,
        )
        cage1 = optimizer.optimize(mol=cage1)

    """

    def __init__(
        self,
        output_dir,
        step_size,
        distance_cut,
        scale_steps=True,
    ):
        """
        Initialize a :class:`Collapser` instance.

        Parameters
        ----------
        output_dir : :class:`str`
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        step_size : :class:`float`
            The relative size of the step to take during collapse.

        distance_cut : :class:`float`
            Distance between distinct building blocks to use as
            threshold for halting collapse in Angstrom.

        scale_steps : :class:`bool`, optional
            Whether to scale the step of each distict building block
            by their relative distance from the molecules centroid.
            Defaults to ``True``

        """

        self._output_dir = output_dir
        self._step_size = step_size
        self._distance_cut = distance_cut
        self._scale_steps = scale_steps

    def _get_inter_bb_distance(self, mol):
        """
        Yield The distances between building blocks in mol.

        Ignores H atoms.

        """

        position_matrix = mol.get_position_matrix()

        for atom1, atom2 in combinations(mol.get_atom_infos(), 2):
            chk1 = (
                atom1.get_atom().get_id() != atom2.get_atom().get_id()
            )
            chk2 = (
                atom1.get_atom().get_atomic_number() != 1
                and atom2.get_atom().get_atomic_number() != 1
            )
            chk3 = (
                atom1.get_building_block_id() !=
                atom2.get_building_block_id()
            )
            if chk1 and chk2 and chk3:
                dist = get_atom_distance(
                    position_matrix=position_matrix,
                    atom1_id=atom1.get_atom().get_id(),
                    atom2_id=atom2.get_atom().get_id()
                )
                yield dist

    def _has_short_contacts(self, mol):
        """
        Calculate if there are short contants in mol.

        """

        return any(
            dist < self._distance_cut
            for dist in self._get_inter_bb_distance(mol)
        )

    def _get_new_position_matrix(self, mol, step, vectors, scales):
        """
        Get the position matrix of the mol after translation.

        """

        new_position_matrix = mol.get_position_matrix()
        for atom in mol.get_atom_infos():
            bb_id = atom.get_building_block_id()
            id = atom.get_atom().get_id()
            pos = mol.get_position_matrix()[id]
            new_position_matrix[id] = (
                pos - step*vectors[bb_id]*scales[bb_id]
            )

        return new_position_matrix

    def _get_bb_vectors(self, mol, bb_atom_ids):
        """
        Get the bb to COM vectors.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        bb_atom_ids : :class:`dict`
            Atom ids (values) in each distinct building block in the
            molecule. Keys are building block ids.

        Returns
        -------
        bb_cent_vectors : :class:`dict`
            Vector from building block (Key is building block id) to
            molecules centroid.

        bb_cent_scales : :class:`dict`
            Relative size of vector between building blocks and
            centroid.

        """

        cent = mol.get_centroid()

        # Get bb COM vector to molecule COM.
        bb_cent_vectors = {
            i: mol.get_centroid(atom_ids=bb_atom_ids[i])-cent
            for i in bb_atom_ids
        }

        # Scale the step size based on the different distances of
        # bbs from the COM. Impacts anisotropic topologies.
        if self._scale_steps:
            max_distance = max(
                np.linalg.norm(bb_cent_vectors[i])
                for i in bb_cent_vectors
            )
            bb_cent_scales = {
                i: np.linalg.norm(bb_cent_vectors[i])/max_distance
                for i in bb_cent_vectors
            }
        else:
            bb_cent_scales = {
                i: 1
                for i in bb_cent_vectors
            }

        return bb_cent_vectors, bb_cent_scales

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`stk.ConstructedMolecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`stk.ConstructedMolecule`
            The optimized molecule.

        """

        # Handle output dir.
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)

        bb_ids = list(set([
            i.get_building_block_id() for i in mol.get_atom_infos()
        ]))
        bb_atom_ids = {i: [] for i in bb_ids}
        for i in mol.get_atom_infos():
            bb_atom_ids[i.get_building_block_id()].append(
                i.get_atom().get_id()
            )

        # Translate each bb along bb_COM_vectors `step`.
        # `step` is the proportion of the bb_COM_vectors that is moved.
        step_no = 0
        step = self._step_size
        while not self._has_short_contacts(mol):
            # Update each step the building block vectors and distance.
            bb_cent_vectors, bb_cent_scales = self._get_bb_vectors(
                mol=mol,
                bb_atom_ids=bb_atom_ids
            )

            new_pos = self._get_new_position_matrix(
                mol=mol,
                step=step,
                vectors=bb_cent_vectors,
                scales=bb_cent_scales
            )
            step_no += 1
            mol = mol.with_position_matrix(new_pos)
            mol.write(
                os.path.join(output_dir, f'collapsed_{step_no}.mol')
            )

        bb_cent_vectors, bb_cent_scales = self._get_bb_vectors(
            mol=mol,
            bb_atom_ids=bb_atom_ids
        )

        # Check that we have not gone too far.
        min_dist = min(
            dist for dist in self._get_inter_bb_distance(mol)
        )
        if min_dist < self._distance_cut / 2:
            # Revert to half the previous step if we have.
            step = -(self._step_size/2)
            new_pos = self._get_new_position_matrix(
                mol=mol,
                step=step,
                vectors=bb_cent_vectors,
                scales=bb_cent_scales
            )
            step_no += 1
            mol = mol.with_position_matrix(new_pos)
            mol.write(
                os.path.join(output_dir, f'collapsed_rev.mol')
            )

        out_file = os.path.join(output_dir, f'collapser.out')
        with open(out_file, 'w') as f:
            f.write(
                f"Collapser algorithm.\n"
                f"====================\n"
                f"Step size: {self._step_size}\n"
                f"Scale steps?: {self._scale_steps}\n"
                f"Distance cut: {self._distance_cut}\n"
                f"====================\n"
                f"Steps run: {step_no}\n"
                f"Minimum inter-bb distance: {min_dist}\n"
                f"====================\n"
            )
        return mol


class Collapserv2(Collapser):
    """
    Collapse stk.ConstructedMolecule to decrease enlarged bonds.

    Smarter optimisation than Collapser.

    """

    def _get_long_bond_infos(self, mol):
        """
        Returns list of long bond infos.

        """

        long_bond_infos = []
        for bond_infos in mol.get_bond_infos():
            if bond_infos.get_building_block() is None:
                long_bond_infos.append(bond_infos)

        return long_bond_infos

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`stk.Molecule`
            The molecule to be optimized.

        Returns
        -------
        mol : :class:`stk.Molecule`
            The optimized molecule.

        """

        # Handle output dir.
        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)

        # Define long bonds to optimise.
        long_bond_infos = self._get_long_bond_infos(mol)
        print(long_bond_infos)
        import sys
        sys.exit()

        BB_ids = list(set([
            i.get_building_block_id() for i in mol.get_atom_infos()
        ]))
        BB_atom_ids = {i: [] for i in BB_ids}
        for i in mol.get_atom_infos():
            BB_atom_ids[i.get_building_block_id()].append(
                i.get_atom().get_id()
            )

        # Translate each BB along BB_COM_vectors `step`.
        # `step` is the proportion of the BB_COM_vectors that is moved.
        step_no = 0
        step = self._step_size
        while not self._has_short_contacts(mol):
            # Update each step the building block vectors and distance.
            BB_cent_vectors, BB_cent_scales = self._get_BB_vectors(
                mol=mol,
                BB_atom_ids=BB_atom_ids
            )

            new_pos = self._get_new_position_matrix(
                mol=mol,
                step=step,
                vectors=BB_cent_vectors,
                scales=BB_cent_scales
            )
            step_no += 1
            mol = mol.with_position_matrix(new_pos)
            mol.write(
                os.path.join(output_dir, f'collapsed_{step_no}.mol')
            )

        BB_cent_vectors, BB_cent_scales = self._get_BB_vectors(
            mol=mol,
            BB_atom_ids=BB_atom_ids
        )

        # Check that we have not gone too far.
        min_dist = min(
            dist for dist in self._get_inter_BB_distance(mol)
        )
        if min_dist < self._distance_cut / 2:
            # Revert to half the previous step if we have.
            step = -(self._step_size/2)
            new_pos = self._get_new_position_matrix(
                mol=mol,
                step=step,
                vectors=BB_cent_vectors,
                scales=BB_cent_scales
            )
            step_no += 1
            mol = mol.with_position_matrix(new_pos)
            mol.write(
                os.path.join(output_dir, f'collapsed_rev.mol')
            )

        out_file = os.path.join(output_dir, f'collapser.out')
        with open(out_file, 'w') as f:
            f.write(
                f"Collapser algorithm.\n"
                f"====================\n"
                f"Step size: {self._step_size}\n"
                f"Scale steps?: {self._scale_steps}\n"
                f"Distance cut: {self._distance_cut}\n"
                f"====================\n"
                f"Steps run: {step_no}\n"
                f"Minimum inter-BB distance: {min_dist}\n"
                f"====================\n"
            )
        return mol
