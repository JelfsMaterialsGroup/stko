import glob
import logging
import os
import sys

import stk
import stko


def main():
    first_line = f"Usage: {__file__}.py"
    if not len(sys.argv) == 2:
        logging.info(f"{first_line} gulp_path")
        sys.exit()
    else:
        gulp_path = sys.argv[1]

    bb1 = stk.BuildingBlock(
        smiles="[Cu+2]",
        functional_groups=(
            stk.SingleAtom(stk.Cu(0, charge=2)) for i in range(4)
        ),
        position_matrix=([0, 0, 0],),
    )
    bb2 = stk.BuildingBlock(
        smiles="O=C(O)c1ccc(Br)cc1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8]~[#1]",
                bonders=(1,),
                deleters=(2,),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#8X1]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )
    structure = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.Paddlewheel(
            metals=bb1,
            ligands=bb2,
            reaction_factory=stk.DativeReactionFactory(
                stk.GenericReactionFactory(
                    bond_orders={
                        frozenset(
                            {
                                stk.GenericFunctionalGroup,
                                stk.SingleAtom,
                            }
                        ): 9,
                    },
                ),
            ),
        ),
    )

    # Perform Gulp optimisation with UFF4MOF.
    # Use conjugate gradient method for a slower, but more stable
    # optimisation.
    gulp_opt = stko.GulpUFFOptimizer(
        gulp_path=gulp_path,
        output_dir="gulp_test_output",
        metal_FF={29: "Cu4+2"},
        conjugate_gradient=True,
        maxcyc=500,
    )
    # Assign the force field.
    gulp_opt.assign_FF(structure)
    # Run optimization.
    structure = gulp_opt.optimize(mol=structure)
    structure.write(os.path.join("gulp_test_output", "opt_structure.mol"))

    target_num_confs = 40
    gulp_MD = stko.GulpUFFMDOptimizer(
        gulp_path=gulp_path,
        metal_FF={29: "Cu4+2"},
        output_dir="gulp_test_output_MD",
        temperature=300,
        equilbration=0.2,
        production=67.9,
        N_conformers=target_num_confs,
        opt_conformers=False,
        save_conformers=True,
    )
    gulp_MD.assign_FF(structure)
    structure = gulp_MD.optimize(structure)
    structure.write(os.path.join("gulp_test_output_MD", "opt_structure.mol"))

    confs_gen = glob.glob("gulp_test_output_MD/conf*.xyz")
    assert len(confs_gen) == target_num_confs


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
