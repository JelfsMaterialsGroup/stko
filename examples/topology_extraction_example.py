import stk
import stko


def main():
    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)C=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage1 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    disconnections = []
    for bi in cage1.get_bond_infos():
        if bi.get_building_block() is None:
            a1id = bi.get_bond().get_atom1().get_id()
            a2id = bi.get_bond().get_atom2().get_id()
            disconnections.append(sorted((a1id, a2id)))

    print(disconnections)
    new_topology_graph = stko.TopologyExtractor()
    new_topology_graph.extract_topology(
        molecule=cage1,
        atom_ids_to_disconnect=disconnections,
    )
    print(new_topology_graph.get_vertex_positions())
    print(new_topology_graph.get_connectivities())


if __name__ == "__main__":
    main()
