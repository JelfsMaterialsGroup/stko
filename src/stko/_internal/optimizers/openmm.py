from typing import Literal

from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper
from openmm import Integrator, LangevinIntegrator, State
from openmm.unit import (
    Quantity,
    angstrom,
    femtoseconds,
    kelvin,
    kilojoule,
    mole,
    nanometer,
    picoseconds,
)

from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT


class OpenMMForceField(Optimizer):
    def __init__(  # noqa: PLR0913
        self,
        force_field: ForceField,
        tolerance: Quantity = 10 * kilojoule / (nanometer * mole),
        max_iterations: int = 0,
        box_vectors: Quantity | None = None,
        define_stereo: bool = False,
        partial_charges_method: Literal["am1bcc", "mmff94"] = "am1bcc",
    ) -> None:
        self._integrator = LangevinIntegrator(
            300 * kelvin, 1 / picoseconds, 0.002 * picoseconds
        )
        self._force_field = force_field
        self._box_vectors = box_vectors
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method
        self._tolerance = tolerance
        self._max_iterations = max_iterations

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
        state: State,
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
