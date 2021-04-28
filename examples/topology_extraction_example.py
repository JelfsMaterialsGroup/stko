import stk
import stko


def main():
    molecule = stk.BuildingBlock.init_from_file(
        'output_directory/poly_uff.mol'
    )

    new_topology_graph = stko.TopologyExtractor()
    print(new_topology_graph)
    new_topology_graph.extract_topology(
        molecule=molecule,
        atom_ids_to_disconnect=((6, 18), ),
        vertex_options={
            1: stk.polymer.vertices._TerminalVertex,
        },
    )
    print([i for i in new_topology_graph.get_vertices()])
    print([i for i in new_topology_graph.get_edges()])


if __name__ == "__main__":
    main()
