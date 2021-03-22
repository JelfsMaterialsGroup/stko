# %%
from dataclasses import dataclass, field
import stk

@dataclass
class Torsion:
    """
    Represents a torsion within a molecule
    
    >>> torsion = Torsion(stk.N(23), stk.C(18), stk.C(19), stk.O(31))
    >>> print(torsion.atom3)
    C(19)
    >>> print(torsion.get_atoms())
    [N(23), C(18), C(19), O(31)]
    """
    atom1: stk.Atom
    atom2: stk.Atom
    atom3: stk.Atom
    atom4: stk.Atom
    
    def get_atoms(self):
        return [self.atom1, self.atom2, self.atom3, self.atom4]
    
    def __iter__(self):
        return iter(self.get_atoms())

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE, verbose=True)

# %%
