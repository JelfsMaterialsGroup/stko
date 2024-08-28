# ruff: noqa: T201
import numpy as np
import stk

import stko


def main() -> None:
    """Run the example."""
    stk_molecule = stk.BuildingBlock("NCCNCCN").with_centroid(
        position=np.array((10, 10, 10))
    )
    rdkit_molecule = stk_molecule.to_rdkit_mol()
    mdanalysis_molecule = stko.MDAnalysis().get_universe(
        mol=stk_molecule,
    )
    print("stk molecule:", stk_molecule)
    print("rdkit molecule:", rdkit_molecule)
    print("Universe:", mdanalysis_molecule)
    print("R_g:", mdanalysis_molecule.atoms.radius_of_gyration())
    print("B_sphere:", mdanalysis_molecule.atoms.bsphere())
    print("Universe COM:", mdanalysis_molecule.atoms.center_of_mass())
    print("stk centroid:", stk_molecule.get_centroid())


if __name__ == "__main__":
    main()
