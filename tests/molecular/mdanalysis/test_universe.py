import stko
import numpy as np


def test_universe(molecule):
    """
    Test :meth:`.MDAnalysis.get_universe`.

    Parameters
    ----------
    molecule : :class:`stk.Molecule`
        The molecule to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = stko.MDAnalysis().get_universe(molecule)

    assert result.atoms.n_atoms == molecule.get_num_atoms()

    for atom, stk_atom in zip(result.atoms, molecule.get_atoms()):
        assert atom.ix == stk_atom.get_id()

    for pos, stk_pos in zip(
        result.atoms.positions, molecule.get_position_matrix()
    ):
        assert np.allclose(pos, stk_pos)
