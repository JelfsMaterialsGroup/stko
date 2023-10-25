"""
Collapser Optimizer
===================

#. :class:`.Collapser`
#. :class:`.CollapserMC`

Optimizer for collapsing enlarged topologies.

"""

import logging
from itertools import combinations
import numpy as np
from collections import defaultdict
from scipy.spatial.distance import pdist
import random
import matplotlib.pyplot as plt
import uuid
import os
import shutil
from stk import PdbWriter

from .optimizers import Optimizer
from ..utilities import get_atom_distance


logger = logging.getLogger(__name__)


class Collapser(Optimizer):
    """
    Collapse stk.ConstructedMolecule to decrease enlarged bonds.

    It is recommended to use the MCHammer version of this code with
    :mod:`MCHammer` [1]_, where a much cleaner version is written.
    The utilities `get_long_bond_ids` will help generate sub units.

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
        optimizer = stko.Collapser(
            output_dir='test_coll',
            step_size=0.05,
            distance_cut=2.0,
            scale_steps=True,
        )
        cage1 = optimizer.optimize(mol=cage1)

    References
    ----------
    .. [1] https://github.com/andrewtarzia/MCHammer

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
            _id = atom.get_atom().get_id()
            pos = mol.get_position_matrix()[_id]
            new_position_matrix[_id] = (
                pos - step*vectors[bb_id]*scales[bb_id]
            )

        return new_position_matrix

    def _get_new_cell_vectors(self, unit_cell, step, scales):
        """
        Get the unit cell vectors after collapse step.

        """

        vector_1 = unit_cell.get_vector_1()
        vector_2 = unit_cell.get_vector_2()
        vector_3 = unit_cell.get_vector_3()

        max_scale = max(list(scales.values()))

        vector_1 = vector_1 - vector_1*step*max_scale
        vector_2 = vector_2 - vector_2*step*max_scale
        vector_3 = vector_3 - vector_3*step*max_scale

        return vector_1, vector_2, vector_3

    def _get_bb_vectors(self, mol, bb_atom_ids):
        """
        Get the building block to COM vectors.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        bb_atom_ids : :class:`dict` mapping :class:`int`: :class:`list`
            Dictionary mapping building block ids (keys) to a list of
            atom ids (values) in each distinct building block in the
            molecule.

        Returns
        -------
        bb_cent_vectors :
            :class:`dict` mapping :class:`int`: :class:`numpy.ndarray`
            Dictionary mapping building block ids (keys) to centroid
            vectors (values) of each distinct building block in the
            molecule.

        bb_cent_scales :
            :class:`dict` mapping :class:`int`: :class:`float`
            Dictionary mapping building block ids (keys) to relative
            magnitude of centroid vectors (values) of each distinct
            building block in the molecule.

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
            norms = {
                i: np.linalg.norm(bb_cent_vectors[i])
                for i in bb_cent_vectors
            }
            max_distance = max(list(norms.values()))
            bb_cent_scales = {
                i: norms[i]/max_distance
                for i in norms
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

        bb_atom_ids = defaultdict(list)
        for i in mol.get_atom_infos():
            bb_atom_ids[i.get_building_block_id()].append(
                i.get_atom().get_id()
            )

        # Translate each building block along bb_COM_vectors by a
        # distance `step`. I.e. `step` is the proportion of the
        # bb_COM_vectors that the building block is moved.
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
                os.path.join(output_dir, 'collapsed_rev.mol')
            )

        out_file = os.path.join(output_dir, 'collapser.out')
        with open(out_file, 'w') as f:
            f.write(
                "Collapser algorithm.\n"
                "====================\n"
                f"Step size: {self._step_size}\n"
                f"Scale steps?: {self._scale_steps}\n"
                f"Distance cut: {self._distance_cut}\n"
                f"====================\n"
                f"Steps run: {step_no}\n"
                f"Minimum inter-bb distance: {min_dist}\n"
                f"====================\n"
            )
        return mol

    def p_optimize(self, mol, unit_cell):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`stk.ConstructedMolecule`
            The molecule to be optimized.

        unit_cell : :class:`.UnitCell`
            The cell to be optimized.

        Returns
        -------
        mol : :class:`stk.ConstructedMolecule`
            The optimized molecule.

        unit_cell : :class:`.UnitCell`
            The optimized cell.

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

        bb_atom_ids = defaultdict(list)
        for i in mol.get_atom_infos():
            bb_atom_ids[i.get_building_block_id()].append(
                i.get_atom().get_id()
            )

        # Translate each building block along bb_COM_vectors by a
        # distance `step`. I.e. `step` is the proportion of the
        # bb_COM_vectors that the building block is moved.
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
                scales=bb_cent_scales,
            )
            step_no += 1
            mol = mol.with_position_matrix(new_pos)
            new_vector_1, new_vector_2, new_vector_3 = (
                self._get_new_cell_vectors(
                    unit_cell=unit_cell,
                    step=step,
                    scales=bb_cent_scales,
                )
            )
            unit_cell = unit_cell.with_cell_from_vectors(
                vector_1=new_vector_1,
                vector_2=new_vector_2,
                vector_3=new_vector_3,
            )
            PdbWriter().write(
                molecule=mol,
                path=os.path.join(
                    output_dir, f'collapsed_{step_no}.pdb'
                ),
                periodic_info=unit_cell,
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
            new_vector_1, new_vector_2, new_vector_3 = (
                self._get_new_cell_vectors(
                    unit_cell=unit_cell,
                    step=step,
                    scales=bb_cent_scales,
                )
            )
            unit_cell = unit_cell.with_cell_from_vectors(
                vector_1=new_vector_1,
                vector_2=new_vector_2,
                vector_3=new_vector_3,
            )
            PdbWriter().write(
                molecule=mol,
                path=os.path.join(
                    output_dir, 'collapsed_rev.pdb'
                ),
                periodic_info=unit_cell,
            )

        out_file = os.path.join(output_dir, 'collapser.out')
        with open(out_file, 'w') as f:
            f.write(
                "Collapser algorithm.\n"
                "====================\n"
                f"Step size: {self._step_size}\n"
                f"Scale steps?: {self._scale_steps}\n"
                f"Distance cut: {self._distance_cut}\n"
                f"====================\n"
                f"Steps run: {step_no}\n"
                f"Minimum inter-bb distance: {min_dist}\n"
                f"====================\n"
            )
        return mol, unit_cell


class CollapserMC(Collapser):
    """
    Collapse molecule to decrease enlarged bonds using MC algorithm.

    It is recommended to use the MCHammer version of this code with
    :mod:`MCHammer` [1]_, where a much cleaner version is written.
    The utilities `get_long_bond_ids` will help generate sub units.

    Smarter optimisation than Collapser using simple Monte Carlo
    algorithm to perform rigid translations of building blocks.

    References
    ----------
    .. [1] https://github.com/andrewtarzia/MCHammer

    """

    def __init__(
        self,
        output_dir,
        step_size,
        target_bond_length,
        num_steps,
        bond_epsilon=50,
        nonbond_epsilon=20,
        nonbond_sigma=1.2,
        nonbond_mu=3,
        beta=2,
        random_seed=None,
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
            The relative size of the step to take during step.

        target_bond_length : :class:`float`
            Target equilibrium bond length for long bonds to minimize
            to.

        num_steps : :class:`int`
            Number of MC moves to perform.

        bond_epsilon : :class:`float`, optional
            Value of epsilon used in the bond potential in MC moves.
            Determines strength of the bond potential.
            Defaults to 50.

        nonbond_epsilon : :class:`float`, optional
            Value of epsilon used in the nonbond potential in MC moves.
            Determines strength of the nonbond potential.
            Defaults to 20.

        nonbond_sigma : :class:`float`, optional
            Value of sigma used in the nonbond potential in MC moves.
            Defaults to 1.2.

        nonbond_mu : :class:`float`, optional
            Value of mu used in the nonbond potential in MC moves.
            Determines the steepness of the nonbond potential.
            Defaults to 3.

        beta : :class:`float`, optional
            Value of beta used in the in MC moves. Beta takes the
            place of the inverse boltzmann temperature.
            Defaults to 2.

        random_seed : :class:`int`, optional
            Random seed to use for MC algorithm. Should only be set
            if exactly reproducible results are required, otherwise
            a system-based random seed should be used for proper
            sampling.

        """

        self._output_dir = output_dir
        self._step_size = step_size
        self._target_bond_length = target_bond_length
        self._num_steps = num_steps
        self._bond_epsilon = bond_epsilon
        self._nonbond_epsilon = nonbond_epsilon
        self._nonbond_sigma = nonbond_sigma
        self._nonbond_mu = nonbond_mu
        self._beta = beta
        if random_seed is None:
            random.seed()
        else:
            random.seed(random_seed)

    def _get_bb_atom_ids(self, mol):

        bb_atom_ids = defaultdict(list)
        for i in mol.get_atom_infos():
            bb_atom_ids[i.get_building_block_id()].append(
                i.get_atom().get_id()
            )

        return bb_atom_ids

    def _get_bond_length(self, mol, bond):

        position_matrix = mol.get_position_matrix()
        return get_atom_distance(
            position_matrix=position_matrix,
            atom1_id=bond.get_atom1().get_id(),
            atom2_id=bond.get_atom2().get_id()
        )

    def _get_bond_vector(self, mol, bond):

        position_matrix = mol.get_position_matrix()
        atom1_pos = position_matrix[bond.get_atom1().get_id()]
        atom2_pos = position_matrix[bond.get_atom2().get_id()]
        return atom2_pos - atom1_pos

    def _get_long_bond_infos(self, mol):
        """
        Returns dict of long bond infos.

        """

        long_bond_infos = {}
        for bond_infos in mol.get_bond_infos():
            if bond_infos.get_building_block() is None:
                ids = (
                    bond_infos.get_bond().get_atom1().get_id(),
                    bond_infos.get_bond().get_atom2().get_id(),
                )
                long_bond_infos[ids] = bond_infos

        return long_bond_infos

    def _get_bb_centroids(self, mol, bb_atom_ids):
        """
        Returns dict of building block centroids.

        """

        bb_centroids = {
            i: mol.get_centroid(atom_ids=bb_atom_ids[i])
            for i in bb_atom_ids
        }

        return bb_centroids

    def _get_cent_to_lb_vector(
        self,
        mol,
        bb_centroids,
        long_bond_infos
    ):
        """
        Returns dict of long bond atom to bb centroid vectors.

        """

        position_matrix = mol.get_position_matrix()
        centroid_to_lb_vectors = {}
        for bb in bb_centroids:
            cent = bb_centroids[bb]
            for b_atom_ids, bond_info in long_bond_infos.items():
                for atom_id in b_atom_ids:
                    atom_info, = mol.get_atom_infos(atom_ids=atom_id)
                    atom_pos = position_matrix[atom_id]
                    if atom_info.get_building_block_id() == bb:
                        centroid_to_lb_vectors[(bb, atom_id)] = (
                            atom_pos - cent,
                        )
                        break

        return centroid_to_lb_vectors

    def _bond_potential(self, distance):
        """
        Define an arbitrary parabolic bond potential.

        This potential has no relation to an empircal forcefield.
        """

        potential = (distance - self._target_bond_length) ** 2
        potential = self._bond_epsilon * potential

        return potential

    def _non_bond_potential(self, distance):
        """
        Define an arbitrary repulsive nonbonded potential.

        This potential has no relation to an empircal forcefield.
        """

        return (
            self._nonbond_epsilon * (
                (self._nonbond_sigma/distance) ** self._nonbond_mu
            )
        )

    def _compute_non_bonded_potential(self, mol):

        # Get all pairwise distances.
        pair_dists = pdist(mol.get_position_matrix())
        nonbonded_potential = np.sum(
            self._non_bond_potential(pair_dists)
        )

        return nonbonded_potential

    def _compute_potential(self, mol, long_bond_infos):

        system_potential = self._compute_non_bonded_potential(mol)
        for long_bond_ids in long_bond_infos:
            long_bond = long_bond_infos[long_bond_ids]

            system_potential += self._bond_potential(
                distance=self._get_bond_length(
                    mol=mol,
                    bond=long_bond.get_bond()
                )
            )

        return system_potential

    def _translate_atoms_along_vector(self, mol, atom_ids, vector):

        new_position_matrix = mol.get_position_matrix()
        for atom in mol.get_atom_infos(atom_ids=atom_ids):
            _id = atom.get_atom().get_id()
            pos = mol.get_position_matrix()[_id]
            new_position_matrix[_id] = pos - vector

        return mol.with_position_matrix(new_position_matrix)

    def _test_move(self, curr_pot, new_pot):

        if new_pot < curr_pot:
            return True
        else:
            exp_term = np.exp(-self._beta*(new_pot-curr_pot))
            rand_number = random.random()

            if exp_term > rand_number:
                return True
            else:
                return False

    def _output_top_lines(self):

        string = (
            '====================================================\n'
            '                Collapser optimisation              \n'
            '                ----------------------              \n'
            '                                                    \n'
            f' step size = {self._step_size} \n'
            f' target bond length = {self._target_bond_length} \n'
            f' num. steps = {self._num_steps} \n'
            f' bond epsilon = {self._bond_epsilon} \n'
            f' nonbond epsilon = {self._nonbond_epsilon} \n'
            f' nonbond sigma = {self._nonbond_sigma} \n'
            f' nonbond mu = {self._nonbond_mu} \n'
            f' beta = {self._beta} \n'
            '====================================================\n\n'
        )

        return string

    def _plot_progess(self, steps, maxds, spots, npots, output_dir):

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(steps, maxds, c='k', lw=2)
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('step', fontsize=16)
        ax.set_ylabel('max long bond length [angstrom]', fontsize=16)
        ax.axhline(y=self._target_bond_length, c='r')
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_dir, 'maxd_vs_step.pdf'),
            dpi=360,
            bbox_inches='tight'
        )
        plt.close()
        # Plot energy vs timestep.
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(steps, spots, c='k', lw=2, label='total')
        ax.plot(steps, npots, c='r', lw=2, label='non bonded')
        # Set number of ticks for x-axis
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.set_xlabel('step', fontsize=16)
        ax.set_ylabel('potential', fontsize=16)
        ax.legend(fontsize=16)
        fig.tight_layout()
        fig.savefig(
            os.path.join(output_dir, 'pot_vs_step.pdf'),
            dpi=360,
            bbox_inches='tight'
        )
        plt.close()

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

        # Define long bonds to optimise.
        long_bond_infos = self._get_long_bond_infos(mol)

        # If no long bonds, then optimisation is done.
        if len(long_bond_infos) == 0:
            return mol

        # Define bb centroid - long bond atom vectors.
        # These are to be maintained during optimisation.
        # centroid_to_lb_vectors = self._get_cent_to_lb_vector(
        #     mol, bb_centroids, long_bond_infos
        # )
        print('###################################################')
        print('WARNING: centroid_to_lb_vectors not maintained yet.')
        print('###################################################')

        system_potential = self._compute_potential(
            mol=mol,
            long_bond_infos=long_bond_infos
        )

        with open(os.path.join(output_dir, 'coll.out'), 'w') as f:
            f.write(self._output_top_lines())
            mol.write(os.path.join(output_dir, 'coll_0.mol'))
            steps = [0]
            passed = []
            spots = [system_potential]
            npots = [self._compute_non_bonded_potential(mol=mol)]
            maxds = [max([
                self._get_bond_length(
                    mol, long_bond_infos[i].get_bond()
                )
                for i in long_bond_infos
            ])]
            f.write(
                'Step system_potential nonbond_potential max_dist '
                'opt_bbs updated?\n'
            )
            f.write(
                f'{steps[-1]} {spots[-1]} {npots[-1]} {maxds[-1]} '
                '-- --\n'
            )
            for step in range(1, self._num_steps):

                # Randomly select a long bond.
                lb_ids = random.choice(list(long_bond_infos.keys()))
                lb_info = long_bond_infos[lb_ids]

                lb_vector = self._get_bond_vector(
                    mol=mol,
                    bond=lb_info.get_bond()
                )

                bb_id_1, bb_id_2 = (
                    i.get_building_block_id()
                    for i in mol.get_atom_infos(atom_ids=lb_ids)
                )

                # Choose bb to move out of the two randomly.
                moving_bb = random.choice([bb_id_1, bb_id_2])
                moving_bb_atom_ids = tuple(
                    i.get_atom().get_id()
                    for i in mol.get_atom_infos()
                    if i.get_building_block_id() == moving_bb
                )

                # Randomly choose between translation along long bond
                # vector or along BB-COM vector.
                # Random number from -1 to 1
                rand = (random.random() - 0.5) * 2

                # Define translation along long bond vector where
                # direction is from force, magnitude is randomly
                # scaled.
                long_bond_translation = (
                    -lb_vector * self._step_size * rand
                )

                # Get bb COM vector to molecule COM.
                cent = mol.get_centroid()
                bb_cent_vector = (
                    mol.get_centroid(atom_ids=moving_bb_atom_ids)-cent
                )
                com_translation = (
                    bb_cent_vector * self._step_size * rand
                )

                translation_vector = random.choice([
                    long_bond_translation,
                    com_translation,
                ])

                # Translate building block.
                # Update atom position of building block.
                mol = self._translate_atoms_along_vector(
                    mol=mol,
                    atom_ids=moving_bb_atom_ids,
                    vector=translation_vector,
                )

                ###################################################
                # Here I want to add rotations
                # To maintain cent_vectors relative orientations.
                # cent_vector_1 = centroid_to_lb_vectors[
                #     (bb_id_1, lb_ids[0])
                # ]
                # cent_vector_2 = centroid_to_lb_vectors[
                #     (bb_id_2, lb_ids[1])
                # ]
                ###################################################

                new_system_potential = self._compute_potential(
                    mol=mol,
                    long_bond_infos=long_bond_infos
                )

                if self._test_move(
                    system_potential,
                    new_system_potential
                ):
                    updated = 'T'
                    system_potential = new_system_potential
                    passed.append(step)
                else:
                    updated = 'F'
                    # Reverse move.
                    mol = self._translate_atoms_along_vector(
                        mol=mol,
                        atom_ids=moving_bb_atom_ids,
                        vector=-translation_vector,
                    )

                mol.write(os.path.join(output_dir, f'coll_{step}.xyz'))
                steps.append(step)
                spots.append(system_potential)
                npots.append(
                    self._compute_non_bonded_potential(mol=mol)
                )
                maxds.append(max([
                    self._get_bond_length(
                        mol, long_bond_infos[i].get_bond()
                    )
                    for i in long_bond_infos
                ]))
                f.write(
                    f'{steps[-1]} {spots[-1]} '
                    f'{npots[-1]} {maxds[-1]} {lb_ids} {updated}\n'
                )
                step += 1

            f.write('\n============================================\n')
            f.write(
                'Optimisation done:\n'
                f'{len(passed)} steps passed: '
                f'{len(passed)/self._num_steps}'
            )

        self._plot_progess(steps, maxds, spots, npots, output_dir)

        return mol
