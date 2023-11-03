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
        print(f"binder distance: {tsa.get_binder_distance(as_building_block)}")
        raise SystemExit("now to do ligand analysis.")

    raise SystemExit("now to do cage analysis : PORE.")
    raise SystemExit("now to do cage analysis : Geom.")
    raise SystemExit("now to do cage analysis : Constructed.")


if __name__ == "__main__":
    main()
