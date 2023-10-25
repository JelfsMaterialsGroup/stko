"""
Molecule Transformer
====================

#. :class:`.MoleculeTransformer`

Class for splitting a molecule into many with new connectors.

"""

import logging
from rdkit.Chem import AllChem as rdkit

import stk

from ..atoms import Du

logger = logging.getLogger(__name__)


class MoleculeTransformer:
    """
    Split an stk.molecule into many with new functional groups.

    Examples
    --------

    Given a molecule, this class allows you to cap a split molecule
    (see :class:`MoleculeSplitter`) at the broken bond with an atom
    defined in `replacer_smarts`.

    .. code-block:: python

        import stk
        import stko

        full_mol = stk.BuildingBlock('C1=CC=NC(=C1)C=NC2=CC=C(C=C2)Br')

        splitter = stko.MoleculeSplitter(
            breaker_smarts='[#6X3]~[#7X2]~[#6X3H1]~[#6X3!H1]',
            bond_deleter_ids=(0, 1),
        )
        split_mols = splitter.split(full_mol)

        transformer = stko.MoleculeTransformer(
            replacer_smarts='[Br]',
            functional_groups=(stk.BromoFactory(), ),
        )
        for split in split_mols:
            transformed_mol = transformer.transform(split)

    """

    def __init__(
        self,
        replacer_smarts,
        functional_groups,
    ):
        """
        Initialize a :class:`.MoleculeTransformer`.

        Parameters
        ----------
        replacer_smarts : :class:`str`
            SMARTS string of atom to replace dummy atoms with. This
            must be a single atom.

        functional_groups : :class:`iterable`
        of :class:`stk.FunctionalGroupFactory`
            Functional group factories to use to define new building
            block.

        Raises
        ------
        :class:`ValueError`
            If `replacer_smarts` does not correspond to a single atom.

        """

        # Add replacers bonded to * atoms.
        _atom_list = list(
            rdkit.MolFromSmarts(replacer_smarts).GetAtoms()
        )
        if len(_atom_list) != 1:
            raise ValueError(
                f'{replacer_smarts} corresponds {len(_atom_list)} '
                'atoms. Should be 1.'
            )

        for i in _atom_list:
            self._replacer = stk.Atom(
                id=0,
                atomic_number=i.GetAtomicNum(),
                charge=0,
            )

        self._functional_groups = functional_groups

    def transform(self, molecule):
        """
        Transform a molecule.

        Parameters
        ----------
        molecule : :class:`stk.Molecule`
            Molecule to modify.

        Returns
        -------
        molecule : :class:`stk.BuildingBlock`
            The resulting molecule.

        """

        rdkit_mol = molecule.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)

        atoms = []
        for a in molecule.get_atoms():
            if isinstance(a, Du):
                atoms.append(stk.Atom(
                    id=a.get_id(),
                    atomic_number=self._replacer.get_atomic_number(),
                    charge=self._replacer.get_charge(),
                ))
            else:
                atoms.append(a)
        atoms = tuple(atoms)

        bonds = tuple(
            stk.Bond(
                atom1=atoms[b.get_atom1().get_id()],
                atom2=atoms[b.get_atom2().get_id()],
                order=b.get_order(),
            )
            for b in molecule.get_bonds()
        )
        position_matrix = molecule.get_position_matrix()
        building_block = stk.BuildingBlock.init(
            atoms=atoms,
            bonds=bonds,
            position_matrix=position_matrix,
        )
        building_block = stk.BuildingBlock.init_from_molecule(
            molecule=building_block,
            functional_groups=self._functional_groups,
        )

        return building_block
