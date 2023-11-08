import os

import numpy as np
import stk
import stko


def main():
    examples_output = "aligner_directory"
    if not os.path.exists(examples_output):
        os.mkdir(examples_output)

    bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles="O=CC(C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            (bb1, bb2),
            optimizer=stk.MCHammer(num_steps=2000),
        ),
    )

    mol_list = [
        (stk.BuildingBlock("NCCNCCN"), (("N", "N"),), True),
        (
            stk.BuildingBlock("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
            (
                ("N", "N"),
                ("O", "O"),
            ),
            True,
        ),
        (stk.BuildingBlock("C1=CN=CN=C1"), (("N", "N"),), True),
        (stk.BuildingBlock("c1ccccc1"), (("C", "C"),), True),
        (stk.BuildingBlock("C1CCCCC1"), (("C", "C"),), True),
        (cage, (("N", "N"),), True),
    ]

    _opt = stko.UFF()
    for i, (mol, pairs, uff_opt) in enumerate(mol_list):
        initial = mol.with_rotation_about_axis(
            1.34,
            np.array((0, 0, 1)),
            np.array((0, 0, 0)),
        )
        mol.write(os.path.join(examples_output, f"unaligned_{i}.mol"))
        initial.write(os.path.join(examples_output, f"init_{i}.mol"))
        aligner = stko.Aligner(
            initial_molecule=initial,
            matching_pairs=pairs,
        )
        if uff_opt:
            mol = aligner.optimize(_opt.optimize(mol))
        else:
            mol = aligner.optimize(mol)
        mol.write(os.path.join(examples_output, f"aligned_{i}.mol"))
        rmsd_calculator = stko.RmsdMappedCalculator(initial)
        rmsd = rmsd_calculator.get_results(mol).get_rmsd()
        print(f"molecule {i} RMSD: {rmsd}")


if __name__ == "__main__":
    main()
