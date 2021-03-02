# %%
from dataclasses import dataclass
from stko.molecular.torsions.torsion import Torsion
import stk
from rdkit.Chem import TorsionFingerprints
from generate_molecule.xor_gate import XorGate

@dataclass
class ConstructedMoleculeTorsioned():
    """
    pass
    """
    stk_molecule: stk.ConstructedMolecule
    
    def __post_init__(self):
        nonring, ring = TorsionFingerprints.CalculateTorsionLists(
            self.stk_molecule.to_rdkit_mol())
        self.torsion_list = [list(atoms[0]) for atoms, ang in nonring]
    
    def get_torsion_list(self):
        'returns a list torsions in the molecule, where each torsion is a list of atom indices'
        return self.torsion_list

    def get_torsions(self):
        'yield the torsions in the molecule'
        for torsion in self.torsion_list:
            yield Torsion(self.stk_molecule.get_atoms(torsion))

    def get_torsions_by_building_block(self):
        'return a set of the torsions corresponding to each building block of this molecule'
        
    def get_torsion_infos(self):
        'yield data about the torsions in the molecule'

    def get_env(self):
        'returns a gym environment corresponding to this molecule'
    
    if __name__ == "__main__":
        xor_gate = ConstructedMoleculeTorsioned(XorGate(3, 8).polymer)
        print(xor_gate.torsion_list)
        

# %%
