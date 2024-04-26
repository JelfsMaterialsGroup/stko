"""Molecular optimisers and property calculators for use with :mod:`stk`."""

from stko import functional_groups, molecule_analysis
from stko._internal.calculators.extractors.orca_extractor import OrcaExtractor
from stko._internal.calculators.open_babel_calculators import OpenBabelEnergy
from stko._internal.calculators.orca_calculators import OrcaEnergy
from stko._internal.calculators.planarity_calculators import (
    PlanarityCalculator,
)
from stko._internal.calculators.rdkit_calculators import MMFFEnergy, UFFEnergy
from stko._internal.calculators.results.energy_results import EnergyResults
from stko._internal.calculators.results.torsion_results import (
    ConstructedMoleculeTorsionResults,
)
from stko._internal.calculators.rmsd_calculators import (
    RmsdCalculator,
    RmsdMappedCalculator,
)
from stko._internal.calculators.shape_calculators import (
    ShapeCalculator,
    ShapeResults,
)
from stko._internal.calculators.torsion_calculators import (
    ConstructedMoleculeTorsionCalculator,
    MatchedTorsionCalculator,
    TorsionCalculator,
)
from stko._internal.calculators.xtb_calculators import XTBEnergy
from stko._internal.molecular.atoms.dummy_atom import Du
from stko._internal.molecular.atoms.positioned_atom import PositionedAtom
from stko._internal.molecular.conversion.md_analysis import MDAnalysis
from stko._internal.molecular.conversion.z_matrix import ZMatrix
from stko._internal.molecular.molecule_modifiers.molecule_splitter import (
    MoleculeSplitter,
)
from stko._internal.molecular.molecule_modifiers.molecule_transformer import (
    MoleculeTransformer,
)
from stko._internal.molecular.networkx.network import Network
from stko._internal.molecular.periodic.unitcell import UnitCell
from stko._internal.molecular.topology_extractor.topology_extractor import (
    TopologyExtractor,
)
from stko._internal.molecular.topology_extractor.topology_info import (
    TopologyInfo,
)
from stko._internal.molecular.torsion.torsion import Torsion
from stko._internal.molecular.torsion.torsion_info import TorsionInfo
from stko._internal.optimizers.aligner import Aligner
from stko._internal.optimizers.collapser import Collapser, CollapserMC
from stko._internal.optimizers.gulp import GulpUFFMDOptimizer, GulpUFFOptimizer
from stko._internal.optimizers.macromodel import (
    MacroModelForceField,
    MacroModelMD,
)
from stko._internal.optimizers.open_babel import OpenBabel
from stko._internal.optimizers.optimizers import (
    Optimizer,
    OptimizerSequence,
    TryCatchOptimizer,
)
from stko._internal.optimizers.rdkit import ETKDG, MMFF, UFF, MetalOptimizer
from stko._internal.optimizers.utilities import MAEExtractor
from stko._internal.optimizers.xtb import XTB, XTBCREST, XTBFF, XTBFFCREST
from stko._internal.types import ConstructedMoleculeT, MoleculeT
from stko._internal.utilities.utilities import get_torsion_info_angles

MoleculeT = MoleculeT  # noqa: PLW0127
"""Type parameter matching any :class:`stk.Molecule` or subclasses."""

ConstructedMoleculeT = ConstructedMoleculeT  # noqa: PLW0127
"""Type parameter matching :class:`stk.ConstructedMolecule` or subclasses."""

__all__ = [
    "functional_groups",
    "molecule_analysis",
    "OrcaExtractor",
    "OpenBabelEnergy",
    "OrcaEnergy",
    "PlanarityCalculator",
    "MMFFEnergy",
    "UFFEnergy",
    "EnergyResults",
    "ConstructedMoleculeTorsionResults",
    "RmsdCalculator",
    "RmsdMappedCalculator",
    "ShapeCalculator",
    "ShapeResults",
    "ConstructedMoleculeTorsionCalculator",
    "MatchedTorsionCalculator",
    "TorsionCalculator",
    "XTBEnergy",
    "Du",
    "PositionedAtom",
    "MDAnalysis",
    "ZMatrix",
    "MoleculeSplitter",
    "MoleculeTransformer",
    "Network",
    "UnitCell",
    "TopologyExtractor",
    "TopologyInfo",
    "Torsion",
    "TorsionInfo",
    "Aligner",
    "Collapser",
    "CollapserMC",
    "GulpUFFMDOptimizer",
    "GulpUFFOptimizer",
    "MacroModelForceField",
    "MacroModelMD",
    "Optimizer",
    "OptimizerSequence",
    "TryCatchOptimizer",
    "ETKDG",
    "MMFF",
    "UFF",
    "MetalOptimizer",
    "MAEExtractor",
    "XTB",
    "XTBFF",
    "XTBCREST",
    "XTBFFCREST",
    "ConstructedMoleculeT",
    "MoleculeT",
    "OpenBabel",
    "get_torsion_info_angles",
]
