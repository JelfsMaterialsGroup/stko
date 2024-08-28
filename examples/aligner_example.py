from pathlib import Path

import numpy as np
import stk

import stko


def main() -> None:
    """Run the example."""
    examples_output = Path("aligner_directory")
    examples_output.mkdir(parents=True, exist_ok=True)

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
        mol.write(examples_output / f"unaligned_{i}.mol")
        initial.write(examples_output / f"init_{i}.mol")
        aligner = stko.Aligner(
            initial_molecule=initial,
            matching_pairs=pairs,
        )
        if uff_opt:
            aligned_mol = aligner.optimize(_opt.optimize(mol))
        else:
            aligned_mol = aligner.optimize(mol)
        aligned_mol.write(examples_output / f"aligned_{i}.mol")
        rmsd_calculator = stko.RmsdMappedCalculator(initial)
        rmsd = rmsd_calculator.get_results(mol).get_rmsd()
        print(f"molecule {i} RMSD: {rmsd}")  # noqa: T201


if __name__ == "__main__":
    main()
