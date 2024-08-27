# ruff: noqa: T201
import logging
from collections import defaultdict

import numpy as np
import stk

import stko


def main() -> None:
    """Run the example."""
    logging.warning(
        "This code is only present in the latest versions of stko that "
        "require Python 3.11!"
    )

    pd = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for _ in range(4)
        ),
        position_matrix=np.array([[0.0, 0.0, 0.0]]),
    )
    ditopic_bb = stk.BuildingBlock(
        smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]", bonders=(1,), deleters=()
            )
        ],
    )
    apdcage = stk.ConstructedMolecule(
        topology_graph=stk.cage.M6L12Cube(
            building_blocks=(pd, ditopic_bb),
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
            optimizer=stk.MCHammer(num_steps=1500),
        ),
    )
    apdcage.write("cage_output/apdcage.mol")
    ligands = stko.molecule_analysis.DecomposeMOC().decompose(
        molecule=apdcage,
        metal_atom_nos=(46,),
    )
    print(f"there are {len(ligands)} ligands")

    tsa = stko.molecule_analysis.DitopicThreeSiteAnalyser()
    ligand_dict = defaultdict(list)
    for id_, lig in enumerate(ligands):
        lig.write(f"cage_output/apdcage_{id_}.mol")
        as_building_block = stk.BuildingBlock.init_from_molecule(
            lig,
            stko.functional_groups.ThreeSiteFactory(smarts="[#6]~[#7X2]~[#6]"),
        )
        ligand_dict["binder_dist"].append(
            tsa.get_binder_distance(as_building_block)
        )
        ligand_dict["binder_binder"].append(
            tsa.get_binder_binder_angle(as_building_block)
        )
        ligand_dict["torsion"].append(
            abs(tsa.get_binder_adjacent_torsion(as_building_block))
        )
        binder_angles = tsa.get_binder_angles(as_building_block)
        ligand_dict["binder_angle"].append(binder_angles[0])
        ligand_dict["binder_angle"].append(binder_angles[1])
        ligand_dict["bite_angle"].append(
            sum(tsa.get_halfbite_angles(as_building_block))
        )

    # Distance between binders in the organic linkers?
    print(f"avg. binder distance: {np.mean(ligand_dict['binder_dist'])}")
    # How twisted is the molecule, or aligned are the binding groups?
    print(
        "avg. binder binder angle: " f"{np.mean(ligand_dict['binder_binder'])}"
    )
    # How twisted is the molecule, or what is the torsion between
    # the binding groups?
    print(
        "avg. binder adjacent torsion: " f"{np.mean(ligand_dict['torsion'])}"
    )
    # What is the angle made by the binders?
    print(f"avg. binder angles: {np.mean(ligand_dict['binder_angle'])}")
    # And the resultant bite-angle [caution!]?
    print(
        "avg. bite angle [caution]: " f"{np.mean(ligand_dict['bite_angle'])}"
    )

    # Get the centroid and atom ids of distinct building blocks.
    # This is similar to the above, but does not perform any disconnections
    # on the constructed molecule or maintain building block chemistry.
    # It simply extracts the parts of building blocks still present in
    # the molecule.
    constructed_analyser = stko.molecule_analysis.ConstructedAnalyser()
    centroids = constructed_analyser.get_building_block_centroids(apdcage)
    print("building block centroids:")
    print(centroids)
    print()

    # Get measures of pore size and cage geometry without any external
    # software. These are just geometrical measures.
    analyser = stko.molecule_analysis.GeometryAnalyser()
    print(
        f"approximate pore size: {analyser.get_min_centroid_distance(apdcage)}"
    )
    print(f"avg cage size: {analyser.get_avg_centroid_distance(apdcage)}")
    m_distances = list(
        analyser.get_metal_distances(
            apdcage,
            metal_atom_nos=(46,),
        ).values()
    )
    print(f"avg. metal distance: {np.mean(m_distances)}")
    m_angles = list(
        analyser.get_metal_centroid_metal_angle(
            apdcage,
            metal_atom_nos=(46,),
        ).values()
    )
    print(f"avg. metal-centroid-angles: {np.mean(m_angles)}")

    # And some geometrical measures.
    print(
        f"avg. N-Pd bond length:"
        f' {np.mean(analyser.calculate_bonds(apdcage)[("N", "Pd")])}'
    )
    print(
        f"N-Pd-N angles: "
        f'{analyser.calculate_angles(apdcage)[("N", "Pd", "N")]}'
    )


if __name__ == "__main__":
    main()
