import os

import stk
import stko


def main():
    examples_output = "splitter_output_directory"
    if not os.path.exists(examples_output):
        os.mkdir(examples_output)

    full_mol = stk.BuildingBlock("C1=CC=NC(=C1)C=NC2=CC=C(C=C2)Br")
    print(full_mol)
    full_mol.write(os.path.join(examples_output, "original.mol"))

    splitter = stko.MoleculeSplitter(
        # This smarts corresponds to an imine bond: `CNCC`.
        breaker_smarts="[#6X3]~[#7X2]~[#6X3H1]~[#6X3!H1]",
        bond_deleter_ids=(0, 1),
    )
    split_mols = splitter.split(full_mol)
    print(split_mols)
    for i, mol in enumerate(split_mols):
        print(mol)
        mol.write(os.path.join(examples_output, f"splits_{i}.xyz"))

    # Test it worked.
    assert len(split_mols) == 2

    # Transform molecule.
    transformer = stko.MoleculeTransformer(
        replacer_smarts="[Br]",
        functional_groups=(stk.BromoFactory(),),
    )
    transformed_mols = tuple(transformer.transform(i) for i in split_mols)

    count = 0
    for i, mol in enumerate(transformed_mols):
        print(mol)
        mol.write(os.path.join(examples_output, f"transform_{i}.mol"))
        count += 1

    # Test it worked.
    assert count == 2

    # Reconstruct.
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(transformed_mols[0], transformed_mols[1]),
            repeating_unit="AB",
            num_repeating_units=1,
        ),
    )
    polymer.write(os.path.join(examples_output, "reconstructed.mol"))

    # Another more complex example.
    full_mol2 = stk.BuildingBlock(
        smiles=(
            "C1=CC=C(C=C1)C#CC2=CC(=C(C=C2C#CC3=CC=CC=C3)C#CC4=CC=CC="
            "C4)C#CC5=CC=CC=C5"
        ),
    )
    print(full_mol2)
    full_mol2.write(os.path.join(examples_output, "original2.mol"))

    splitter = stko.MoleculeSplitter(
        breaker_smarts="[#6]~[#6]#[#6]",
        bond_deleter_ids=(0, 1),
    )
    split_mols2 = splitter.split(full_mol2)

    for i, mol in enumerate(split_mols2):
        print(mol)
        mol.write(os.path.join(examples_output, f"splits2_{i}.xyz"))

    # Transform molecule.
    transformer = stko.MoleculeTransformer(
        replacer_smarts="[Br]",
        functional_groups=(stk.BromoFactory(),),
    )
    transformed_mols2 = tuple(transformer.transform(i) for i in split_mols2)

    for i, mol in enumerate(transformed_mols2):
        print(mol)
        mol.write(os.path.join(examples_output, f"transform2_{i}.mol"))


if __name__ == "__main__":
    main()
