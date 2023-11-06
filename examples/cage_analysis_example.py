import stk
import stko


def main():
    pd = stk.BuildingBlock(
        smiles="[Pd+2]",
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
        ),
        position_matrix=[[0.0, 0.0, 0.0]],
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
    print(len(ligands))

    tsa = stko.molecule_analysis.DitopicThreeSiteAnalyser()
    for i, lig in enumerate(ligands):
        lig.write(f"cage_output/apdcage_{i}.mol")
        as_building_block = stk.BuildingBlock.init_from_molecule(
            lig,
            stko.functional_groups.ThreeSiteFactory(smarts="[#6]~[#7X2]~[#6]"),
        )

        # Distance between binders in the organic linkers?
        print(f"binder distance: {tsa.get_binder_distance(as_building_block)}")

        # How twisted is the molecule, or aligned are the binding groups?
        print(
            "binder binder angle: "
            f"{tsa.get_binder_binder_angle(as_building_block)}"
        )

        # How twisted is the molecule, or what is the torsion between
        # the binding groups?
        print(
            "binder adjacent torsion: "
            f"{abs(tsa.get_binder_adjacent_torsion(as_building_block))}"
        )

        # What is the angle made by the binders?
        print(f"binder angles: {tsa.get_binder_angles(as_building_block)}")

        # And the resultant bite-angle [caution!]?
        print(
            "bite angle [caution]: "
            f"{sum(tsa.get_halfbite_angles(as_building_block))}"
        )

        raise SystemExit("now to do ligand analysis.")

    # Get the centroid and atom ids of distinct building blocks.
    # This is similar to the above, but does not perform any disconnections
    # on the constructed molecule or maintain building block chemistry.
    # It simply extracts the parts of building blocks still present in
    # the molecule.
    analyser = stko.molecule_analysis.ConstructedAnalyser()
    atom_ids = analyser.get_building_block_atom_ids(apdcage)
    print(atom_ids)
    centroids = analyser.get_building_block_centroids(apdcage)
    print(centroids)

    raise SystemExit("now to do cage analysis : PORE.")
    raise SystemExit("now to do cage analysis : Geom.")


if __name__ == "__main__":
    main()
