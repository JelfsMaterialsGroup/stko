# ruff: noqa: T201
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

    # Zipping calculations together.
    mols = [
        # Perfect linear short.
        stk.BuildingBlock("C#C"),
        # Kind of linear long.
        stk.BuildingBlock("C#CC#C"),
        stk.BuildingBlock("CCCC"),
        # Flat.
        stk.BuildingBlock("C1=CC=CC=C1"),
        # Less flat.
        stk.BuildingBlock("C1CCCCC1"),
        # Sphere - CH4.
        stk.BuildingBlock("C"),
        # Sphere - adamantane.
        stk.BuildingBlock("C1C3CC2CC(CC1C2)C3"),
        # Constructed Molecule.
        polymer,
    ]
    shape_calc = stko.ShapeCalculator()
    for mol in mols:
        print(
            mol,
            stko.ShapeResults(
                shape_calc.calculate(mol)
            ).get_spherocity_index(),
            stko.ShapeResults(shape_calc.calculate(mol)).get_eccentricity(),
        )


if __name__ == "__main__":
    main()
