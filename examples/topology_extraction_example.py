import stk
import stko


def main():
    molecule = stk.BuildingBlock.init_from_file(
        'output_directory/poly_uff.mol'
    )

    new_topology_graph = stko.TopologyExtractor(
        vertex_options={
            2: stk.polymer.vertices._LinearVertex,
            1: stk.polymer.vertices._TerminalVertex,
        },
    )
    print(new_topology_graph)
    new_topology_graph.extract_topology(
        molecule=molecule,
        bond_ids_to_break=((6, 18), ),
    )
    print([i for i in new_topology_graph.get_vertices()])
    print([i for i in new_topology_graph.get_edges()])
    print([i for i in new_topology_graph.get_disconnected_graphs()])


if __name__ == "__main__":
    main()
