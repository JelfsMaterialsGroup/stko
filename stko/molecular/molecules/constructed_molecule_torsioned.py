# %%
from dataclasses import dataclass
from collections import defaultdict
from stko.molecular.torsions.torsion_info import TorsionInfo
from stko.molecular.torsions.torsion import Torsion
import stk
from rdkit.Chem import TorsionFingerprints

class ConstructedMoleculeTorsioned():
    """
    
    """
    stk_molecule: stk.ConstructedMolecule
    
    def __init__(self, stk_molecule: stk.ConstructedMolecule):
        self.stk_molecule = stk_molecule.clone()
        nonring, ring = TorsionFingerprints.CalculateTorsionLists(
            self.stk_molecule.to_rdkit_mol())
        self.torsions = [Torsion(*atoms[0]) for atoms, ang in nonring]
        self.set_atom_maps()
    
    def get_torsion_list(self):
        'returns a list torsions in the molecule, where each torsion is a list of atom indices'
        return [torsion.get_atoms() for torsion in self.torsions]

    def get_torsions(self):
        'yield the torsions in the molecule'
        yield from self.torsions

    def get_torsion_infos_by_building_block(self):
        'return a set of the torsions corresponding to each building block of this molecule'
        torsion_infos_by_building_block = defaultdict(list)
        for torsion_info in self.get_torsion_infos():
            torsion_infos_by_building_block[torsion_info.building_block_id].append(torsion_info)
        return torsion_infos_by_building_block
        
    def get_torsion_infos(self):
        'yield data about the torsions in the molecule'
        for torsion, atom_ids in zip(self.get_torsions(), self.get_torsion_list()):
            atom_infos = list(self.stk_molecule.get_atom_infos(atom_ids))
            building_block_ids = {atom_info.get_building_block_id() for atom_info in atom_infos}
            if len(building_block_ids) == 1:
                bb_atoms = [atom_info.get_building_block_atom() for atom_info in atom_infos]
                building_block_torsion = Torsion(*bb_atoms)
                building_block = atom_infos[0].get_building_block()
                building_block_id = building_block_ids.pop()
                yield TorsionInfo(torsion, building_block_torsion, building_block, building_block_id)
            else:
                yield TorsionInfo(torsion, None, None, None)
    
    def set_torsions(self, torsions):
        'sets the torsions'
        self.torsions = torsions
        
    def transfer_torsions(self, building_block_map):
        """
        building_block_map is a map from each building block of this molecule
        to a ConstructedMoleculeTorsioned which contains its torsions
        """
        self.torsions = []
        print('Here!!!!!')
        print(f'self.get_building_blocks() = {self.get_building_blocks()}')
        for id, building_block in self.get_building_blocks().items():
            print(f'id: {id}\nbuilding_block: {building_block}')
            for building_block_torsion in building_block_map[building_block].get_torsions():
                torsion = Torsion([self.atom_maps[id][atom] for atom in building_block_torsion])
                self.torsions.append(torsion)
                
    def get_building_blocks(self):
        return {atom_info.get_building_block_id(): atom_info.get_building_block()
                for atom_info in self.stk_molecule.get_atom_infos()}
     
    def set_atom_maps(self):
        """
        map from building block atom ids to constructed molecule atoms for a
        specified building block id
        """
        self.atom_maps = defaultdict(dict)
        for atom_info in self.stk_molecule.get_atom_infos():
            current_atom_map = self.atom_maps[atom_info.get_building_block_id()]
            current_atom_map[atom_info.get_building_block_atom()] = atom_info.get_atom()
        return self.atom_maps
    
if __name__ == "__main__":
    from generate_molecule.xor_gate import XorGate
    xor_gate = ConstructedMoleculeTorsioned(XorGate(3, 8).polymer)
    print(xor_gate.get_torsion_list())
    print([torsion_info.building_block_id for torsion_info in xor_gate.get_torsion_infos()])
    # print([(key, value.building_block_torsion) 
    #       for key, value in xor_gate.get_torsion_infos_by_building_block().items()])
    print(xor_gate.get_torsion_infos_by_building_block())
        

# %%
