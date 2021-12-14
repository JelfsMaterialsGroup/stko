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

from ..atoms import Du

logger = logging.getLogger(__name__)


class MoleculeSplitter:
    """
    Split an stk.molecule into many with dummy atoms.

    Examples
    --------

    Given a molecule, this class allows you to break bonds based on
    `breaker_smarts` between the atoms in `bond_deleter_ids`.

    .. code-block:: python

        import stk
        import stko

        full_mol = stk.BuildingBlock('C1=CC=NC(=C1)C=NC2=CC=C(C=C2)Br')

        splitter = stko.MoleculeSplitter(
            breaker_smarts='[#6X3]~[#7X2]~[#6X3H1]~[#6X3!H1]',
            bond_deleter_ids=(0, 1),
        )
        split_mols = splitter.split(full_mol)

    """

    def __init__(
        self,
        breaker_smarts,
        bond_deleter_ids,
    ):
        """
        Initialize a :class:`.MoleculeSplitter`.

        Parameters
        ----------
        breaker_smarts : :class:`str`
            SMARTS string used to find the substructure to break.

        bond_deleter_ids : :class:`tuple` of :class:`int`
            Index of atoms in `breaker_smarts` to break bond between.

        """

        self._breaker_smarts = breaker_smarts
        self._bond_deleter_ids = bond_deleter_ids

    def split(self, molecule):
        """
        Split a molecule.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            Molecule to modify.

        Returns
        -------
        molecules : :class:`iterable` of :class:`stk.BuildingBlock`
            The resulting list of molecules.

        """

        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        bond_pair_ids_to_delete = []
        reactable_atom_ids = []
        for atom_ids in rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmarts(self._breaker_smarts)
        ):
            # Get the atom ids of the query on either end of the broken
            # bond, translated from query ids to molecule atom ids.
            translated_bond_atom_ids = (
                atom_ids[self._bond_deleter_ids[0]],
                atom_ids[self._bond_deleter_ids[1]],
            )
            # Get tuples of the bonds to break by atom ids.
            bond_pair_ids_to_delete.extend(
                tuple(sorted((i, j)))
                for i, j in combinations(translated_bond_atom_ids, 2)
            )
            # Get the atom ids of the atoms on either end of the broken
            # bond at which reactions can occur.
            reactable_atom_ids.extend(
                i for i in translated_bond_atom_ids
            )

        # Get the rdkit molecule bond ids associated with the bonds to
        # delete.
        bond_ids_to_delete = []
        for bond in rdkit_mol.GetBonds():
            idxs = tuple(sorted((
                bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            )))
            if idxs in bond_pair_ids_to_delete:
                bond_ids_to_delete.append(bond.GetIdx())

        # Delete bonds and atoms.
        fragments = rdkit.GetMolFrags(
            rdkit.FragmentOnBonds(
                mol=rdkit_mol,
                bondIndices=tuple(bond_ids_to_delete),
                addDummies=True,
            ),
            asMols=True,
        )

        molecules = []
        for frag in fragments:
            atoms = []
            for a in frag.GetAtoms():
                if a.GetAtomicNum() == 0:
                    atoms.append(Du(a.GetIdx()))
                else:
                    atoms.append(
                        stk.Atom(
                            id=a.GetIdx(),
                            atomic_number=a.GetAtomicNum(),
                            charge=a.GetFormalCharge(),
                        )
                    )
            atoms = tuple(atoms)

            bonds = tuple(
                stk.Bond(
                    atom1=atoms[b.GetBeginAtomIdx()],
                    atom2=atoms[b.GetEndAtomIdx()],
                    order=(
                        9 if b.GetBondType() == rdkit.BondType.DATIVE
                        else b.GetBondTypeAsDouble()
                    )
                )
                for b in frag.GetBonds()
            )
            position_matrix = frag.GetConformer().GetPositions()
            molecules.append(
                stk.BuildingBlock.init(
                    atoms=atoms,
                    bonds=bonds,
                    position_matrix=position_matrix,
                )
            )

        return molecules
