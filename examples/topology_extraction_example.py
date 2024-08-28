# ruff: noqa: T201
from pathlib import Path

import stk

import stko


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles="O=CC(C=O)C=O",
        functional_groups=[stk.AldehydeFactory()],
    )
    cage1 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    broken_bonds_by_id = []
    disconnectors: set[int] = set()
    for bi in cage1.get_bond_infos():
        if bi.get_building_block() is None:
            a1id = bi.get_bond().get_atom1().get_id()
            a2id = bi.get_bond().get_atom2().get_id()
            a, b = sorted((a1id, a2id))
            broken_bonds_by_id.append((a, b))
            disconnectors.update((a1id, a2id))

    print(broken_bonds_by_id)
    print(disconnectors)
    print("--")
    new_topology_graph = stko.TopologyExtractor()
    tg_info = new_topology_graph.extract_topology(
        molecule=cage1,
        broken_bonds_by_id=broken_bonds_by_id,
        disconnectors=disconnectors,
    )
    print(tg_info.get_vertex_positions())
    print(tg_info.get_connectivities())
    print(tg_info.get_edge_pairs())
    output_directory = Path("output_directory")
    cage1.write(output_directory / "tg_cage.mol")
    tg_info.write(output_directory / "tg_info.pdb")


if __name__ == "__main__":
    main()
