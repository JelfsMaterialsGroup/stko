# ruff: noqa: T201
import sys
from pathlib import Path

import stk

import stko


def main() -> None:
    """Run the example."""
    if len(sys.argv) > 1:
        orca_path = sys.argv[1]
    else:
        print("usage: orca_example.py orca_path")
        sys.exit()

    bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=(0, 0),
            num_repeating_units=1,
        )
    )

    examples_output = Path("orca_output_directory")
    examples_output.mkdir(parents=True, exist_ok=True)

    # Run optimisations.
    etkdg = stko.ETKDG()
    polymer = etkdg.optimize(polymer)

    orca_ey_1 = stko.OrcaEnergy(
        orca_path=orca_path,
        topline="! SP B97-3c",
        basename="example1",
        output_dir=f"{examples_output}/orca_e1_dir",
    )
    print(orca_ey_1.get_energy(polymer))

    uff = stko.UFF()
    polymer = uff.optimize(polymer)
    orca_ey_2 = stko.OrcaEnergy(
        orca_path=orca_path,
        topline="! SP B97-3c",
        basename="example2",
        output_dir=f"{examples_output}/orca_e2_dir",
    )
    print(orca_ey_2.get_energy(polymer))

    orca_ey_3 = stko.OrcaEnergy(
        orca_path=orca_path,
        topline="! SP B97-3c",
        basename="example3",
        output_dir=f"{examples_output}/orca_e3_dir",
        write_input_only=True,
    )
    orca_ey_3.get_results(polymer)


if __name__ == "__main__":
    main()
