# ruff: noqa: T201

import logging
import time
from pathlib import Path

import openmm
import stk
from openff.toolkit import ForceField

import stko


def integrator(
    *,
    temperature: float,
    friction: float,
    time_step: float,
) -> openmm.LangevinIntegrator:
    """Define integrator."""
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    integrator.setRandomNumberSeed(127)
    return integrator


def main() -> None:
    """Run the example."""
    output_directory = Path("cage_openmm_example_directory")
    output_directory.mkdir(exist_ok=True)

    # Construct a cage.
    bb1 = stk.BuildingBlock(
        smiles="C1CCC(C(C1)N)N",
        functional_groups=[stk.PrimaryAminoFactory()],
    )
    bb2 = stk.BuildingBlock(
        smiles="C1=C(C=C(C=C1C=O)C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )
    cage.write(output_directory / "unopt_cage.mol")

    # Settings.
    force_field = ForceField("openff_unconstrained-2.1.0.offxml")
    partial_charges = "espaloma-am1bcc"
    temperature = 700 * openmm.unit.kelvin
    friction = 10 / openmm.unit.picoseconds
    time_step = 1 * openmm.unit.femtoseconds

    # Define sequence.
    optimisation_sequence = stko.OptimizerSequence(
        # Restricted true to optimised the constructed bonds.
        stko.OpenMMForceField(
            force_field=force_field,
            restricted=True,
            partial_charges_method=str(partial_charges),
        ),
        # Unrestricted optimisation.
        stko.OpenMMForceField(
            # Load the openff-2.1.0 force field appropriate for
            # vacuum calculations (without constraints)
            force_field=force_field,
            restricted=False,
            partial_charges_method=partial_charges,
        ),
        # Molecular dynamics, short for equilibration.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output_directory / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method=partial_charges,
            # Frequency here is not related to the num confs tested.
            reporting_freq=100,
            trajectory_freq=100,
            # 10 ps
            num_steps=10_000,
            num_conformers=10,
            platform="CUDA",
            conformer_optimiser=stko.OpenMMForceField(
                force_field=force_field,
                restricted=False,
                partial_charges_method=partial_charges,
            ),
        ),
        # Long MD, for collecting lowest energy conformers.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output_directory / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method=partial_charges,
            # Frequency here is not related to the num confs tested.
            reporting_freq=100,
            trajectory_freq=100,
            # 0.2 ns
            num_steps=200_000,
            # 1 every 4 ps
            num_conformers=50,
            platform="CUDA",
            conformer_optimiser=stko.OpenMMForceField(
                force_field=force_field,
                restricted=False,
                partial_charges_method=partial_charges,
            ),
        ),
    )
    logging.info("starting calculations...")
    st = time.time()
    optimised_cage = optimisation_sequence.optimize(cage)
    et = time.time()
    logging.info("calculations done!")

    optimised_cage.write(output_directory / "opt_cage.mol")
    logging.info(
        "cage energy: %s kJmol-1 in %s s",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="espaloma-am1bcc",
        ).get_energy(optimised_cage),
        round(et - st, 2),
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
