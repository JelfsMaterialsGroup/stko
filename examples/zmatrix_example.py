import os

import stk
import stko


def main():
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

    examples_output = "output_directory"
    if not os.path.exists(examples_output):
        os.mkdir(examples_output)

    print(stko.ZMatrix().get_zmatrix(bb1))
    bb1.write(os.path.join(examples_output, "bb1.mol"))
    print(stko.ZMatrix().get_zmatrix(bb2))
    bb2.write(os.path.join(examples_output, "bb2.mol"))
    print(stko.ZMatrix().get_zmatrix(polymer))
    polymer.write(os.path.join(examples_output, "polymer.mol"))


if __name__ == "__main__":
    main()
