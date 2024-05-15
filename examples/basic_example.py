from pathlib import Path

import stk
import stko


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCNCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=(0, 0),
            num_repeating_units=1,
        )
    )

    examples_output = Path("output_directory")
    if not examples_output.exists():
        examples_output.mkdir()

    # Run optimisations.
    uff = stko.UFF()
    polymer = uff.optimize(polymer)
    polymer.write(examples_output / "poly_uff.mol")
    mmff = stko.MMFF()
    polymer = mmff.optimize(polymer)
    polymer.write(examples_output / "poly_mmff.mol")
    etkdg = stko.ETKDG()
    polymer = etkdg.optimize(polymer)
    polymer.write(examples_output / "poly_etkdg.mol")


if __name__ == "__main__":
    main()
