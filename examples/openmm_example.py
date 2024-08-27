# ruff: noqa: T201

from pathlib import Path

import openmm
import stk
from openff.toolkit import ForceField

import stko


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("C1CCC(C(C1)N)N", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        "C1=C(C=C(C=C1C=O)C=O)C=O", [stk.AldehydeFactory()]
    )
    output = Path("openmm_example_directory")
    output.mkdir(exist_ok=True, parents=True)
    cage = stk.ConstructedMolecule(stk.cage.FourPlusSix([bb1, bb2]))
    cage.write(output / "unopt_cage.mol")
    ff_optimizer = stko.OpenMMForceField(
        # Load the openff-2.1.0 force field appropriate for
        # vacuum calculations (without constraints)
        force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
        partial_charges_method="mmff94",
        max_iterations=10,
    )
    ff_cage = ff_optimizer.optimize(cage)
    ff_cage.write(output / "ff_opt_cage.mol")
    print(
        "ff_opt_cage",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="mmff94",
        ).get_energy(ff_cage),
    )

    ff_optimizer = stko.OpenMMForceField(
        # Load the openff-2.1.0 force field appropriate for
        # vacuum calculations (without constraints)
        force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
        restricted=False,
        partial_charges_method="mmff94",
    )
    ff_cage = ff_optimizer.optimize(cage)
    ff_cage.write(output / "ffunrest_opt_cage.mol")
    print(
        "ffunrest_opt_cage",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="mmff94",
        ).get_energy(ff_cage),
    )

    ff_optimizer = stko.OpenMMForceField(
        # Load the openff-2.1.0 force field appropriate for
        # vacuum calculations (without constraints)
        force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
        restricted=True,
        partial_charges_method="mmff94",
    )
    ff_cage = ff_optimizer.optimize(cage)
    ff_cage.write(output / "ffrest_opt_cage.mol")
    print(
        "ffrest_opt_cage",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="mmff94",
        ).get_energy(ff_cage),
    )

    temperature = 750 * openmm.unit.kelvin
    friction = 10 / openmm.unit.picoseconds
    time_step = 1 * openmm.unit.femtoseconds
    integrator = openmm.LangevinIntegrator(temperature, friction, time_step)
    integrator.setRandomNumberSeed(127)
    md_optimizer = stko.OpenMMMD(
        force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
        output_directory=output / "md_optimisation",
        integrator=integrator,
        random_seed=275,
        partial_charges_method="mmff94",
        reporting_freq=5,
        trajectory_freq=5,
        num_steps=10000,
        num_conformers=50,
        platform="CUDA",
    )
    md_cage = md_optimizer.optimize(cage)
    md_cage.write(output / "md_opt_cage.mol")
    print(
        "md_opt_cage",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="mmff94",
        ).get_energy(md_cage),
    )


if __name__ == "__main__":
    main()
