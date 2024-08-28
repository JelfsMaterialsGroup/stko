import pathlib
import shutil
from copy import copy
from typing import Literal, Protocol

import rdkit.Chem as rdkit  # noqa: N813
import stk
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper
from openmm import app, openmm

from stko._internal.calculators.openmm_calculators import OpenMMEnergy
from stko._internal.optimizers.optimizers import NullOptimizer, Optimizer
from stko._internal.types import MoleculeT
from stko._internal.utilities.exceptions import InputError
from stko._internal.utilities.utilities import get_atom_distance


class EnergyCalculator(Protocol):
    def get_energy(self, mol: stk.Molecule) -> float: ...


class OpenMMForceField(Optimizer):
    """Uses OpenMM to optimise molecules.

    .. tip::

        You can get all force fields with:

        .. testcode::

            import openff.toolkit
            openff.toolkit.typing.engines.smirnoff.get_available_force_fields()

    Parameters:
        force_field:
            The force field to use.
        restricted:
            If ``True`` then an optimization is performed only on bonds
            created during the `ConstructedMolecule`
            creation.
            All building block bonds will be fixed.
            If ``False`` then all bonds are optimized.
        tolerance:
            The energy tolerance to which the system should be minimized
        max_iterations:
            The maximum number of iterations to perform. If this is 0,
            minimization is continued until the results converge without
            regard to how many iterations it takes.
        box_vectors:
            The unit-wrapped box vectors of this topology.
        define_stereo:
            Toggle calculation of stereochemistry.
        partial_charges_method:
            The method to use for calculating partial charges.
            The default ``"am1bcc"`` is semi-empirical and may be slow.

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
        partial_charges_method: Literal[
            "am1bcc", "mmff94", "gasteiger", "am1-mulliken"
        ] = "am1bcc",
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
        # Handle issue with existing context.
        integrator = copy(self._integrator)

        rdkit_mol = mol.to_rdkit_mol()
        if self._define_stereo:
            rdkit.AssignStereochemistry(rdkit_mol)

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
        simulation = app.Simulation(topology, system, integrator)
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
    """Optimise a molecule with OpenMM and Molecular Dynamics.

    Parameters:

        force_field:
            The force field to use.

        output_directory:
            The directory to which the output files should be written.

        reporting_freq:
            How often the simulation properties should be written in
            time steps.

        trajectory_freq:
            How often the trajectory should be written in time steps.

        integrator:
            The integrator to use, :class:`openmm.openmm.LangevinIntegrator`
            by default.

        num_steps:
            The number of steps to simulate.

        num_conformers:
            The number of conformers to sample during the MD run.

        box_vectors:
            The box vectors to use.

        define_stereo:
            Toggle calculation of stereochemistry.

        partial_charges_method:
            The method to use for calculating partial charges.
            The default ``"am1bcc"`` is semi-empirical and may be slow.

        random_seed:
            The random seed to use.

        initial_temperature:
            The initial temperature to use.

        platform:
            The platform to use.

        conformer_optimiser:
            The optimiser to use for the conformers.

        energy_calculator:
            The energy calculator to use to evaluate conformers.

    """

    def __init__(  # noqa: PLR0913
        self,
        force_field: ForceField,
        output_directory: pathlib.Path,
        reporting_freq: int,
        trajectory_freq: int,
        integrator: openmm.Integrator | None = None,
        num_steps: int = 100,
        num_conformers: int = 50,
        box_vectors: openmm.unit.Quantity | None = None,
        define_stereo: bool = False,
        partial_charges_method: Literal[
            "am1bcc", "mmff94", "gasteiger", "am1-mulliken"
        ] = "am1bcc",
        random_seed: int = 108,
        initial_temperature: openmm.unit.Quantity = 300 * openmm.unit.kelvin,
        platform: Literal["CUDA"] | None = None,
        conformer_optimiser: Optimizer | None = None,
        energy_calculator: EnergyCalculator | None = None,
    ) -> None:
        self._output_directory = output_directory
        self._trajectory_data = self._output_directory / "trajectory_data.dat"
        self._trajectory_file = (
            self._output_directory / "trajectory_structures.pdb"
        )

        if integrator is None:
            integrator = openmm.LangevinIntegrator(
                300 * openmm.unit.kelvin,
                1 / openmm.unit.picoseconds,
                0.25 * openmm.unit.femtoseconds,
            )
            integrator.setRandomNumberSeed(34)

        self._integrator = integrator
        self._random_seed = random_seed
        self._num_steps = num_steps
        self._num_conformers = num_conformers

        self._force_field = force_field
        self._box_vectors = box_vectors
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method
        self._initial_temperature = initial_temperature
        self._reporting_freq = reporting_freq
        self._trajectory_freq = trajectory_freq

        if conformer_optimiser is None:
            conformer_optimiser = NullOptimizer()

        self._conformer_optimiser = conformer_optimiser

        if energy_calculator is None:
            energy_calculator = OpenMMEnergy(
                force_field=force_field,
                box_vectors=box_vectors,
                define_stereo=define_stereo,
                partial_charges_method=partial_charges_method,
            )
        self._energy_calculator = energy_calculator

        if platform is not None:
            self._platform = openmm.Platform.getPlatformByName(platform)
            self._properties: dict[str, str] | None = {
                "CudaPrecision": "mixed"
            }
        else:
            self._platform = None
            self._properties = None

    def _add_trajectory_reporter(
        self,
        simulation: app.Simulation,
    ) -> app.Simulation:
        simulation.reporters.append(
            app.PDBReporter(
                file=str(self._trajectory_file),
                reportInterval=self._trajectory_freq,
            )
        )
        return simulation

    def _add_reporter(self, simulation: app.Simulation) -> app.Simulation:
        simulation.reporters.append(
            app.StateDataReporter(
                file=str(self._trajectory_data),
                reportInterval=self._reporting_freq,
                step=True,
                potentialEnergy=True,
                kineticEnergy=True,
                totalEnergy=False,
                temperature=True,
                volume=False,
                density=False,
                progress=False,
                remainingTime=False,
                speed=False,
                totalSteps=self._num_steps,
                separator=",",
            )
        )
        return simulation

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        if self._output_directory.exists():
            shutil.rmtree(self._output_directory)
        self._output_directory.mkdir(parents=True)

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
        simulation = interchange.to_openmm_simulation(
            integrator=self._integrator,
            platform=self._platform,
            platformProperties=self._properties,
        )
        simulation.context.setVelocitiesToTemperature(
            self._initial_temperature,
            self._random_seed,
        )
        simulation.minimizeEnergy()

        # Add reporters.
        simulation = self._add_reporter(simulation=simulation)
        simulation = self._add_trajectory_reporter(simulation=simulation)

        chunk_size = self._num_steps // self._num_conformers
        min_energy = float("inf")
        min_energy_conformer = None
        for _ in range(self._num_conformers):
            simulation.step(chunk_size)
            state = simulation.context.getState(
                getPositions=True, getEnergy=True
            )
            conformer_mol = self._update_stk_molecule(mol, state)

            conformer_mol = self._conformer_optimiser.optimize(conformer_mol)

            energy = self._energy_calculator.get_energy(conformer_mol)

            if energy < min_energy:
                min_energy = energy
                min_energy_conformer = conformer_mol

        return min_energy_conformer  # type: ignore[return-value]

    def _update_stk_molecule(
        self,
        molecule: MoleculeT,
        state: openmm.State,
    ) -> MoleculeT:
        positions = state.getPositions(asNumpy=True)
        return molecule.with_position_matrix(positions * 10)
