# ruff: noqa: S101
import argparse
import logging
from pathlib import Path

import numpy as np
import stk

import stko


def main() -> None:
    """Run the example."""
    args = _parse_args()

    iron_atom = stk.BuildingBlock(
        smiles="[Fe+2]",
        functional_groups=(
            stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(6)
        ),
        position_matrix=np.array([[0, 0, 0]]),
    )
    bb2 = stk.BuildingBlock(
        smiles="C1=NC(C=NBr)=CC=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )
    mcomplex = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.OctahedralDelta(
            metals=iron_atom,
            ligands=bb2,
            optimizer=stk.MCHammer(),
        ),
    )
    # Assign Bromo functional groups to the metal complex.
    iron_oct_delta = stk.BuildingBlock.init_from_molecule(
        molecule=mcomplex,
        functional_groups=[stk.BromoFactory()],
    )

    # Define building blocks.
    bb3 = stk.BuildingBlock(
        smiles=("C1=CC(=CC=C1C2=CC=C(C=C2)Br)Br"),
        functional_groups=[stk.BromoFactory()],
    )

    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M4L6TetrahedronSpacer(
            building_blocks={
                iron_oct_delta: (0, 1, 2, 3),
                bb3: (4, 5, 6, 7, 8, 9),
            },
            optimizer=stk.MCHammer(),
        ),
    )

    # Perform Gulp optimisation with UFF4MOF.
    # Use conjugate gradient method for a slower, but more stable
    # optimisation.
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=args.gulp_path,
        output_dir="gulp_test_output",
        metal_FF={26: "Fe4+2"},
        conjugate_gradient=True,
        maxcyc=500,
    )
    # Assign the force field.
    gulp_opt.assign_FF(cage)
    # Run optimization.
    structure = gulp_opt.optimize(mol=cage)
    structure.write(Path("gulp_test_output") / "opt_structure.mol")

    target_num_confs = 40
    gulp_md = stko.GulpUFFMDOptimizer(
        gulp_path=args.gulp_path,
        metal_FF={26: "Fe4+2"},
        output_dir="gulp_test_output_MD",
        temperature=300,
        equilbration=0.2,
        production=67.9,
        N_conformers=target_num_confs,
        opt_conformers=False,
        save_conformers=True,
    )
    gulp_md.assign_FF(structure)
    structure = gulp_md.optimize(structure)
    structure.write(Path("gulp_test_output_MD") / "opt_structure.mol")

    confs_gen = list(Path("gulp_test_output_MD").glob("conf*.xyz"))
    assert len(confs_gen) == target_num_confs


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("gulp_path", type=str)
    return parser.parse_args()


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
