# This is an example for incorporating the benzene fixture
import stk
import stko
from fixtures.benzene import benzene_build

def test_benzene(benzene_build):
    mol = benzene_build
    optimizer = stko.ETKDG()
    opt_benzene = optimizer.optimize(mol)
    pos = [x for x in opt_benzene.get_atomic_positions()]
    bond_length = ((pos[0] - pos[1]) ** 2).sum() ** 0.5
    assert abs(1.4 - bond_length) < 0.05