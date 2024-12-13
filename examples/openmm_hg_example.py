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
    output = Path("openmm_hg_example_directory")
    output.mkdir(exist_ok=True, parents=True)

    input_cage = stk.BuildingBlock.init_from_file(
        Path("cage_openmm_example_directory/opt_cage.mol")
    )

    guest1 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("CC"),
        displacement=(0.0, 3.0, 0.0),
    )
    guest2 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("C1CCCC1"),
    )
    hg_complex = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(input_cage),
            guests=(guest1, guest2),
            optimizer=stk.Spinner(),
        ),
    )
    hg_complex.write(output / "unopt_hgcomplex.mol")

    # Settings.
    force_field = ForceField("openff_unconstrained-2.1.0.offxml")
    temperature = 300 * openmm.unit.kelvin
    friction = 10 / openmm.unit.picoseconds
    time_step = 1 * openmm.unit.femtoseconds

    # Define sequence.
    optimisation_sequence = stko.OptimizerSequence(
        # Unrestricted optimisation.
        stko.OpenMMForceField(
            # Load the openff-2.1.0 force field appropriate for
            # vacuum calculations (without constraints)
            force_field=force_field,
            restricted=False,
            partial_charges_method="espaloma-am1bcc",
        ),
        # Molecular dynamics, short for equilibration.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method="espaloma-am1bcc",
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
                partial_charges_method="espaloma-am1bcc",
            ),
        ),
        # Long MD, for collecting lowest energy conformers.
        stko.OpenMMMD(
            force_field=force_field,
            output_directory=output / "md_optimisation",
            integrator=integrator(
                temperature=temperature,
                friction=friction,
                time_step=time_step,
            ),
            random_seed=275,
            partial_charges_method="espaloma-am1bcc",
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
                partial_charges_method="espaloma-am1bcc",
            ),
        ),
    )
    logging.info("starting calculations...")
    st = time.time()
    optimised_complex = optimisation_sequence.optimize(hg_complex)
    et = time.time()
    logging.info("calculations done!")

    optimised_complex.write(output / "opt_complex.mol")
    logging.info(
        "hgcomplex energy: %s kJmol-1 in %s s",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="espaloma-am1bcc",
        ).get_energy(optimised_complex),
        round(et - st, 2),
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
