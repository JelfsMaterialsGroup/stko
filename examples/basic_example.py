import stk
import stko

import os
import shutil


def main():
    bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=[0, 0],
            num_repeating_units=1
        )
    )

    examples_output = 'output_directory'
    if os.path.exists(examples_output):
        shutil.rmtree(examples_output)

    os.mkdir(examples_output)

    # Run optimisations.
    uff = stko.UFF()
    polymer = uff.optimize(polymer)
    polymer.write(f'{examples_output}/poly_uff.mol')
    mmff = stko.MMFF()
    polymer = mmff.optimize(polymer)
    polymer.write(f'{examples_output}/poly_mmff.mol')
    etkdg = stko.ETKDG()
    polymer = etkdg.optimize(polymer)
    polymer.write(f'{examples_output}/poly_etkdg.mol')



if __name__ == "__main__":
    main()