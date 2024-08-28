# ruff: noqa: T201
from pathlib import Path

import numpy as np
import stk

import stko


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock("O=CCC=O", [stk.AldehydeFactory()])
    polymer = stk.ConstructedMolecule(
        stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit="AB",
            orientations=(0, 0),
            num_repeating_units=1,
        )
    )
    bb2 = stk.BuildingBlock("O=CC(C=O)C=O", [stk.AldehydeFactory()])
    cage = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )
    cage2 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(bb1, bb2),
            optimizer=stk.MCHammer(),
        ),
    )

    # Produce a Fe+2 atom with 6 functional groups.
    iron_atom = stk.BuildingBlock(
        smiles="[Fe+2]",
        functional_groups=(
            stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(6)
        ),
        position_matrix=np.array([[0, 0, 0]]),
    )

    # Define coordinating ligand with dummy bromine groups and
    # metal coordinating functional groups.
    oct_bb = stk.BuildingBlock(
        smiles="C1=NC(C=NBr)=CC=C1",
        functional_groups=[
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#35]",
                bonders=(1,),
                deleters=(),
            ),
            stk.SmartsFunctionalGroupFactory(
                smarts="[#6]~[#7X2]~[#6]",
                bonders=(1,),
                deleters=(),
            ),
        ],
    )

    # Build iron complex with delta stereochemistry.
    iron_oct_delta_bb = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.OctahedralDelta(
            metals=iron_atom,
            ligands=oct_bb,
            optimizer=stk.MCHammer(),
        ),
    )

    # Assign Bromo functional groups to the metal complex.
    iron_oct_delta = stk.BuildingBlock.init_from_molecule(
        molecule=iron_oct_delta_bb,
        functional_groups=[stk.BromoFactory()],
    )

    # Define spacer building block.
    bb3 = stk.BuildingBlock(
        smiles=("C1=CC(C2=CC=C(Br)C=C2)=C" "C=C1Br"),
        functional_groups=[stk.BromoFactory()],
    )

    # Build an M4L6 Tetrahedron with a spacer.
    moc = stk.ConstructedMolecule(
        topology_graph=stk.cage.M4L6TetrahedronSpacer(
            building_blocks=(
                iron_oct_delta,
                bb3,
            ),
            optimizer=stk.MCHammer(),
        ),
    )

    host = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="NC1CCCCC1N",
                    functional_groups=[
                        stk.PrimaryAminoFactory(),
                    ],
                ),
                stk.BuildingBlock(
                    smiles="O=Cc1cc(C=O)cc(C=O)c1",
                    functional_groups=[stk.AldehydeFactory()],
                ),
            ),
            optimizer=stk.MCHammer(),
        ),
    )
    guest1 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("BrBr"),
        displacement=(0.0, 3.0, 0.0),
    )
    guest2 = stk.host_guest.Guest(
        building_block=stk.BuildingBlock("C1CCCC1"),
    )

    hg_complex = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=stk.BuildingBlock.init_from_molecule(host),
            guests=(guest1, guest2),
        ),
    )

    cycle = stk.ConstructedMolecule(
        topology_graph=stk.macrocycle.Macrocycle(
            building_blocks=(
                stk.BuildingBlock(
                    smiles="[Br]CC[Br]",
                    functional_groups=[stk.BromoFactory()],
                ),
            ),
            repeating_unit="A",
            num_repeating_units=8,
            optimizer=stk.MCHammer(),
        ),
    )
    axle = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(
                stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
                stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()]),
            ),
            repeating_unit="AB",
            num_repeating_units=3,
            optimizer=stk.MCHammer(),
        )
    )
    rotaxane = stk.ConstructedMolecule(
        topology_graph=stk.rotaxane.NRotaxane(
            axle=stk.BuildingBlock.init_from_molecule(axle),
            cycles=(stk.BuildingBlock.init_from_molecule(cycle),),
            repeating_unit="A",
            num_repeating_units=1,
        ),
    )

    examples_output = Path("output_directory")
    examples_output.mkdir(parents=True, exist_ok=True)

    structures = [
        ("bb1", bb1),
        ("bb2", bb2),
        ("bb3", bb3),
        ("polymer", polymer),
        ("cage", cage),
        ("cage2", cage2),
        ("m_iron_oct_delta", iron_oct_delta),
        ("m_moc", moc),
        ("hg_complex", hg_complex),
        ("cycle", cycle),
        ("axle", axle),
        ("rotaxane", rotaxane),
    ]

    # Run optimisations.
    for name, struct in structures:
        if "m_" in name:
            rdkit_uff_energy = None
        else:
            rdkit_uff_energy = stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(struct)
        obabel_uff_energy = stko.OpenBabelEnergy("uff").get_energy(struct)
        struct.write(examples_output / f"obabel_{name}_unopt.mol")
        opt = stko.OpenBabel(
            forcefield="uff",
            repeat_steps=10,
            sd_steps=20,
            cg_steps=20,
        )
        opt_struct = opt.optimize(struct)
        opt_struct.write(examples_output / f"obabel_{name}_opt.mol")

        if "m_" in name:
            new_rdkit_uff_energy = None
        else:
            new_rdkit_uff_energy = stko.UFFEnergy(
                ignore_inter_interactions=False
            ).get_energy(opt_struct)
        new_obabel_uff_energy = stko.OpenBabelEnergy("uff").get_energy(
            opt_struct
        )
        print(
            f"{name}:\n"
            f"rdkit: {rdkit_uff_energy}, obabel: {obabel_uff_energy}\n"
            f"opt rdkit: {new_rdkit_uff_energy}, "
            f"opt obabel: {new_obabel_uff_energy}\n"
        )


if __name__ == "__main__":
    main()
