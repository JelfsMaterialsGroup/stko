import stk
import numpy as np
import stko


def main():

    an_stk_molecule = stk.BuildingBlock('NCCNCCN').with_centroid(
        position=np.array((10, 10, 10))
    )
    an_rdkit_molecule = an_stk_molecule.to_rdkit_mol()
    an_mdanalysis_molecule = stko.MDAnalysis().get_universe(
        mol=an_stk_molecule,
    )
    print('stk molecule:', an_stk_molecule)
    print('rdkit molecule:', an_rdkit_molecule)
    print('Universe:', an_mdanalysis_molecule)
    print('R_g:', an_mdanalysis_molecule.atoms.radius_of_gyration())
    print('B_sphere:', an_mdanalysis_molecule.atoms.bsphere())
    print(
        'Universe COM:',
        an_mdanalysis_molecule.atoms.center_of_mass()
    )
    print('stk centroid:', an_stk_molecule.get_centroid())


if __name__ == "__main__":
    main()
