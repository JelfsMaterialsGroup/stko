import stk
import numpy as np
import stko

import os


def main():
    examples_output = 'output_directory'
    if not os.path.exists(examples_output):
        os.mkdir(examples_output)

    mol = stk.BuildingBlock('NCCNCCN')
    mol2 = mol.with_rotation_about_axis(
        1.34, np.array((0, 0, 1)), np.array((0, 0, 0)),
    )
    mol.write(os.path.join(examples_output, 'unaligned.mol'))
    mol2.write(os.path.join(examples_output, 'init.mol'))
    aligner = stko.Aligner(
        initial_molecule=mol2,
        matching_pairs=(('N', 'N'), ),
    )
    mol = aligner.optimize(mol)
    mol.write(os.path.join(examples_output, 'aligned.mol'))


if __name__ == "__main__":
    main()
