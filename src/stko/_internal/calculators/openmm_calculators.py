import logging
from collections import abc
from copy import copy
from typing import Literal

import stk
from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule, RDKitToolkitWrapper
from openmm import app, openmm

logger = logging.getLogger(__name__)


class OpenMMEnergy:
    """Uses OpenMM to calculate energy.

    Parameters:

        force_field:
            The force field to use.

        box_vectors:
            The box vectors to use.

        define_stereo:
            Toggle calculation of stereochemistry.

        partial_charges_method:
            The method to use for calculating partial charges.
            The default ``"am1bcc"`` is semi-empirical and may be slow.

    """

    def __init__(
        self,
        force_field: ForceField,
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

        self._force_field = force_field
        self._box_vectors = box_vectors
        self._define_stereo = define_stereo
        self._partial_charges_method = partial_charges_method

    def calculate(self, mol: stk.Molecule) -> abc.Generator:
        # Handle issue with existing context.
        integrator = copy(self._integrator)

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

        # Define simulation.
        simulation = app.Simulation(topology, system, integrator)
        # Set positions from structure.
        simulation.context.setPositions(
            mol.get_position_matrix() * openmm.unit.angstrom
        )
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )

        yield state.getPotentialEnergy().in_units_of(
            openmm.unit.kilojoules_per_mole
        )

    def get_energy(self, mol: stk.Molecule) -> float:
        """Calculate the energy of `mol` in kilojoules per mole.

        Parameters:
            mol:
                The molecule for which to calculate energy.

        Returns:
            The energy.

        """
        return next(self.calculate(mol)).value_in_unit(
            openmm.unit.kilojoules_per_mole
        )
