from typing import Literal

from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper
from openmm import Integrator, LangevinIntegrator, State
from openmm.app import Simulation
from openmm.unit import Quantity, angstrom, kelvin, picosecond, picoseconds

from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT


class OpenMMForceField(Optimizer):
    def __init__(  # noqa: PLR0913
        self,
        force_field: ForceField,
        box_vectors: Quantity | None = None,
        integrator: Integrator | None = None,
        num_steps: int = 10_000,
        define_stereo: bool = False,
        partial_charges_method: Literal["am1bcc", "mmff94"] = "am1bcc",
    ) -> None:
        if integrator is None:
            integrator = LangevinIntegrator(
                300 * kelvin, 1 / picosecond, 0.002 * picoseconds
            )
        self._integrator = integrator
        self._force_field = force_field
        self._box_vectors = box_vectors
        self._num_steps = num_steps
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method

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
