from typing import Literal

import stk
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper
from openmm import app, openmm

from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT
from stko._internal.utilities.exceptions import InputError
from stko._internal.utilities.utilities import get_atom_distance


class OpenMMForceField(Optimizer):
    """Uses OpenMM to optimise molecules.

    Parameters:

        restricted:
            If ``True`` then an optimization is performed only on bonds
            created during the `ConstructedMolecule`
            creation.
            All building block bonds will be fixed.
            If ``False`` then all bonds are optimized.

    """

    def __init__(  # noqa: PLR0913
        self,
        force_field: ForceField,
        restricted: bool = False,
        tolerance: openmm.unit.Quantity = 10
        * openmm.unit.kilojoule
        / (openmm.unit.nanometer * openmm.unit.mole),
        max_iterations: int = 0,
        box_vectors: openmm.unit.Quantity | None = None,
        define_stereo: bool = False,
        partial_charges_method: Literal["am1bcc", "mmff94"] = "am1bcc",
        platform: Literal["CUDA"] | None = None,
    ) -> None:
        self._integrator = openmm.LangevinIntegrator(
            300 * openmm.unit.kelvin,
            1 / openmm.unit.picoseconds,
            0.002 * openmm.unit.picoseconds,
        )
        self._restricted = restricted
        self._force_field = force_field
        self._box_vectors = box_vectors
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method
        self._tolerance = tolerance
        self._max_iterations = max_iterations

        if platform is not None:
            self._platform = openmm.Platform.getPlatformByName(platform)
            if platform == "CUDA":
                self._properties = {"CudaPrecision": "mixed"}
            else:
                self._properties = None
        else:
            self._platform = None
            self._properties = None

    def _add_atom_constraints(
        self,
        system: openmm.System,
        molecule: stk.Molecule,
    ) -> openmm.System:
        if not isinstance(molecule, stk.ConstructedMolecule):
            msg = (
                f"{molecule} needs to be a ConstructedMolecule for "
                f"`restricted`==`True` (it is {self._restricted})"
            )
            raise InputError(msg)

        pos_mat = molecule.get_position_matrix()
        for bond_info in molecule.get_bond_infos():
            # Intra bb bond, so skip.
            if bond_info.get_building_block_id() is None:
                continue
            bond = bond_info.get_bond()
            atom1_id = bond.get_atom1().get_id()
            atom2_id = bond.get_atom2().get_id()

            system.addConstraint(
                particle1=atom1_id,
                particle2=atom2_id,
                distance=get_atom_distance(
                    position_matrix=pos_mat,
                    atom1_id=atom1_id,
                    atom2_id=atom2_id,
                )
                * openmm.unit.angstrom,
            )

        return system

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        rdkit_mol = mol.to_rdkit_mol()
        if self._define_stereo:
            pass

        molecule = Molecule.from_rdkit(
            rdmol=rdkit_mol,
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
        )

        if self._partial_charges_method == "mmff94":
            molecule.assign_partial_charges(
                self._partial_charges_method,
                toolkit_registry=RDKitToolkitWrapper(),
            )

        topology = molecule.to_topology()
        if self._box_vectors is not None:
            topology.box_vectors = self._box_vectors

        interchange = Interchange.from_smirnoff(
            force_field=self._force_field,
            topology=topology,
            positions=mol.get_position_matrix() * openmm.unit.angstrom,
            charge_from_molecules=[molecule],
        )
        system = interchange.to_openmm_system()
        # Add constraints.
        if self._restricted:
            system = self._add_atom_constraints(system, mol)

        # Define simulation.
        simulation = app.Simulation(
            topology,
            system,
            self._integrator,
            platform=self._platform,
            platformProperties=self._properties,
        )
        # Set positions from structure.
        simulation.context.setPositions(
            mol.get_position_matrix() * openmm.unit.angstrom
        )

        simulation.minimizeEnergy(
            tolerance=self._tolerance,
            maxIterations=self._max_iterations,
        )
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )
        return self._update_stk_molecule(mol, state)

    def _update_stk_molecule(
        self,
        molecule: MoleculeT,
        state: openmm.State,
    ) -> MoleculeT:
        positions = state.getPositions(asNumpy=True)
        return molecule.with_position_matrix(positions * 10)


class OpenMMMD(Optimizer):
    def __init__(  # noqa: PLR0913
        self,
        force_field: ForceField,
        integrator: Integrator | None = None,
        num_steps: int = 100,
        box_vectors: Quantity | None = None,
        define_stereo: bool = False,
        partial_charges_method: Literal["am1bcc", "mmff94"] = "am1bcc",
        random_seed: int = 108,
        initial_temperature: Quantity = 300 * kelvin,
    ) -> None:
        if integrator is None:
            integrator = LangevinIntegrator(
                300 * kelvin, 1 / picoseconds, 0.25 * femtoseconds
            )
            integrator.setRandomNumberSeed(34)
        self._integrator = integrator
        self._num_steps = num_steps
        self._force_field = force_field
        self._box_vectors = box_vectors
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method
        self._random_seed = random_seed
        self._initial_temperature = initial_temperature

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        rdkit_mol = mol.to_rdkit_mol()
        if self._define_stereo:
            pass

        molecule = Molecule.from_rdkit(
            rdmol=rdkit_mol,
            allow_undefined_stereo=True,
            hydrogens_are_explicit=True,
        )

        if self._partial_charges_method == "mmff94":
            molecule.assign_partial_charges(
                self._partial_charges_method,
                toolkit_registry=RDKitToolkitWrapper(),
            )

        topology = molecule.to_topology()
        if self._box_vectors is not None:
            topology.box_vectors = self._box_vectors

        interchange = Interchange.from_smirnoff(
            force_field=self._force_field,
            topology=topology,
            positions=mol.get_position_matrix() * angstrom,
            charge_from_molecules=[molecule],
        )
        simulation = interchange.to_openmm_simulation(self._integrator)
        simulation.context.setVelocitiesToTemperature(
            self._initial_temperature, self._random_seed
        )
        simulation.minimizeEnergy()
        simulation.step(self._num_steps)
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )
        return self._update_stk_molecule(mol, state)

    def _update_stk_molecule(
        self,
        molecule: MoleculeT,
        state: State,
    ) -> MoleculeT:
        positions = state.getPositions(asNumpy=True)
        return molecule.with_position_matrix(positions * 10)
