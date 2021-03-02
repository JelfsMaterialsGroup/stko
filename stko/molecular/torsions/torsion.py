# %%
from dataclasses import dataclass
import stk

@dataclass
class Torsion:
    """
    Represents a torsion within a molecule
    
    >>> torsion = Torsion(stk.N(23), stk.C(18), stk.C(19), stk.O(31))
    >>> print(torsion.atom1)
    N(23)
    >>> print(torsion.atom3)
    C(19)
    """
    atom1: stk.Atom
    atom2: stk.Atom
    atom3: stk.Atom
    atom4: stk.Atom

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE, verbose=True)

# %%
