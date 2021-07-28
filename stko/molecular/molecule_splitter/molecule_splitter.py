"""
Molecule Splitter
=================

#. :class:`.MoleculeSplitter`

Class for splitting a molecule into many with new connectors.

"""

import logging
from rdkit.Chem import AllChem as rdkit
from itertools import combinations

import stk

logger = logging.getLogger(__name__)


class MoleculeSplitter:
    """
    Split an stk.molecule into many with new functional groups.

    Examples
    --------

    Given a molecule, this class allows you to break bonds based on
    `breaker_smarts` between the atoms in `bond_deleter_ids` and cap
    the broken bond by an atom defined in `replacer_smarts`.

    .. code-block:: python

        import stk
        import stko

        full_mol = stk.BuildingBlock('C1=CC=NC(=C1)C=NC2=CC=C(C=C2)Br')

        splitter = stko.MoleculeSplitter(
            breaker_smarts='[#6X3]~[#7X2]~[#6X3H1]~[#6X3!H1]',
            bond_deleter_ids=(0, 1),
            replacer_smarts='[Br]',
            functional_groups=(stk.BromoFactory(), ),
        )
        split_mols = splitter.split(full_mol)

    """

    def __init__(
        self,
        breaker_smarts,
        bond_deleter_ids,
        replacer_smarts,
        functional_groups,
    ):
        """
        Initialize a :class:`MoleculeSplitter`.

        Parameters
        ----------
        breaker_smarts : :class:`str`
            SMARTS string to find the substructure to break.

        bond_deleter_ids : :class:`tuple` of :class:`int`
            Index of atoms in `breaker_smarts` to breka bond between.

        replacer_smarts : :class:`str`
            SMARTS string of atom to replace dummy atoms with. This
            must be a single atom.

        functional_groups :
        :class:`iterable` of :class:`stk.FunctionalGroupFactory`
            Functional group factories to use to define new building
            block.

        """

        self._breaker_smarts = breaker_smarts
        self._bond_deleter_ids = bond_deleter_ids
        self._replacer_smarts = replacer_smarts
        self._functional_groups = functional_groups

    def split(self, molecule):
        """
        Split a molecule into multiple molecules.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            Molecule to split.

        Returns
        -------
        :class:`iterable` of :class:`stk.BuildingBlock`
            The resulting list of building blocks, with new functional
            groups assigned.

        """

        new_building_blocks = []

        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        bond_pair_ids_to_delete = []
        reactable_atom_ids = []
        for atom_ids in rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmarts(self._breaker_smarts)
        ):
            translated_bond_atom_ids = (
                atom_ids[self._bond_deleter_ids[0]],
                atom_ids[self._bond_deleter_ids[1]],
            )
            bond_pair_ids_to_delete.extend(
                tuple(sorted((i, j)))
                for i, j in combinations(translated_bond_atom_ids, 2)
            )
            reactable_atom_ids.extend(
                i for i in translated_bond_atom_ids
            )

        bond_ids_to_delete = []
        for bond in rdkit_mol.GetBonds():
            idxs = tuple(sorted((
                bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            )))
            if idxs in bond_pair_ids_to_delete:
                bond_ids_to_delete.append(bond.GetIdx())

        # Delete bonds and atoms.
        fragments = rdkit.GetMolFrags(rdkit.FragmentOnBonds(
            mol=rdkit_mol,
            bondIndices=tuple(bond_ids_to_delete),
            addDummies=True,
        ), asMols=True)

        # Add replacers bonded to * atoms.
        for i in rdkit.MolFromSmarts(self._replacer_smarts).GetAtoms():
            replacer_atom = rdkit.Atom(i.GetAtomicNum())
            replacer_atom.SetFormalCharge(0)
            break

        for frag in fragments:
            editable = rdkit.EditableMol(frag)
            for atom in frag.GetAtoms():
                # 0 if dummy atom.
                if atom.GetAtomicNum() == 0:
                    editable.ReplaceAtom(atom.GetIdx(), replacer_atom)

            frag_rdkit = editable.GetMol()
            rdkit.Kekulize(frag_rdkit)
            new_building_blocks.append(
                stk.BuildingBlock.init_from_rdkit_mol(
                    molecule=frag_rdkit,
                    functional_groups=self._functional_groups
                )
            )

        return new_building_blocks
