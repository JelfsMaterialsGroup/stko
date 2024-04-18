"""Molecular optimisers and property calculators for use with ``stk``."""

from stko import functional_groups, molecule_analysis
from stko._internal.calculators.extractors.orca_extractor import OrcaExtractor
from stko._internal.calculators.extractors.xtb_extractor import XTBExtractor
from stko._internal.calculators.open_babel_calculators import OpenBabelEnergy
from stko._internal.calculators.orca_calculators import OrcaEnergy
from stko._internal.calculators.planarity_calculators import (
    PlanarityCalculator,
)
from stko._internal.calculators.rdkit_calculators import MMFFEnergy, UFFEnergy
from stko._internal.calculators.results.energy_results import EnergyResults
from stko._internal.calculators.results.orca_results import OrcaResults
from stko._internal.calculators.results.planarity_results import (
    PlanarityResults,
)
from stko._internal.calculators.results.rmsd_results import RmsdResults
from stko._internal.calculators.results.shape_results import ShapeResults
from stko._internal.calculators.results.torsion_results import (
    ConstructedMoleculeTorsionResults,
    TorsionResults,
)
from stko._internal.calculators.results.xtb_results import XTBResults
from stko._internal.calculators.rmsd_calculators import (
    RmsdCalculator,
    RmsdMappedCalculator,
)
from stko._internal.calculators.shape_calculators import ShapeCalculator
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
from stko._internal.molecular.periodic.utilities import (
    cap_absolute_value,
    get_approximate_cell_size,
    get_from_parameters,
)
from stko._internal.molecular.topology_extractor.topology_extractor import (
    TopologyExtractor,
)
from stko._internal.molecular.topology_extractor.topology_info import (
    TopologyInfo,
)
from stko._internal.molecular.torsion.torsion import Torsion
from stko._internal.molecular.torsion.torsion_info import TorsionInfo
from stko._internal.optimizers.aligner import (
    Aligner,
    AlignmentPotential,
)
from stko._internal.optimizers.collapser import (
    Collapser,
    CollapserMC,
)
from stko._internal.optimizers.gulp import (
    GulpUFFMDOptimizer,
    GulpUFFOptimizer,
)
from stko._internal.optimizers.macromodel import (
    MacroModel,
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
from stko._internal.utilities.exceptions import (
    CalculatorError,
    ConvergenceError,
    ConversionError,
    DifferentAtomError,
    DifferentMoleculeError,
    ExpectedMetalError,
    ForceFieldError,
    ForceFieldSetupError,
    InputError,
    InvalidSolventError,
    LewisStructureError,
    NotCompletedError,
    NotStartedError,
    OptimizerError,
    PathError,
    SettingConflictError,
    WrapperNotInstalledError,
)
from stko._internal.utilities.utilities import (
    calculate_angle,
    calculate_dihedral,
    get_atom_distance,
    get_torsion_info_angles,
    is_valid_xtb_solvent,
    unit_vector,
    vector_angle,
)

__all__ = [
    "functional_groups",
    "molecule_analysis",
    "is_valid_xtb_solvent",
    "get_atom_distance",
    "get_torsion_info_angles",
    "calculate_angle",
    "calculate_dihedral",
    "vector_angle",
    "unit_vector",
    "Aligner",
    "AlignmentPotential",
    "Collapser",
    "CollapserMC",
    "ExpectedMetal",
    "GulpUFFMDOptimizer",
    "GulpUFFOptimizer",
    "MacroModel",
    "MacroModelForceField",
    "MacroModelMD",
    "OpenBabel",
    "Optimizer",
    "OptimizerSequence",
    "TryCatchOptimizer",
    "MMFF",
    "UFF",
    "ETKDG",
    "MetalOptimizer",
    "XTB",
    "XTBCREST",
    "XTBFF",
    "XTBFFCREST",
    "MAEExtractor",
    "PositionedAtom",
    "Du",
    "MDAnalysis",
    "ZMatrix",
    "MoleculeSplitter",
    "MoleculeTransformer",
    "Network",
    "UnitCell",
    "cap_absolute_value",
    "get_approximate_cell_size",
    "get_from_parameters",
    "TopologyExtractor",
    "TopologyInfo",
    "Torsion",
    "TorsionInfo",
    "XTBEnergy",
    "TorsionCalculator",
    "ConstructedMoleculeTorsionCalculator",
    "MatchedTorsionCalculator",
    "ShapeCalculator",
    "RmsdCalculator",
    "RmsdMappedCalculator",
    "MMFFEnergy",
    "UFFEnergy",
    "OrcaEnergy",
    "OpenBabelEnergy",
    "PlanarityCalculator",
    "XTBResults",
    "TorsionResults",
    "ConstructedMoleculeTorsionResults",
    "ShapeResults",
    "RmsdResults",
    "PlanarityResults",
    "OrcaResults",
    "EnergyResults",
    "XTBExtractor",
    "OrcaExtractor",
    "WrapperNotInstalledError",
    "OptimizerError",
    "ForceFieldError",
    "ForceFieldSetupError",
    "CalculatorError",
    "DifferentAtomError",
    "DifferentMoleculeError",
    "ConvergenceError",
    "ConversionError",
    "ExpectedMetalError",
    "PathError",
    "LewisStructureError",
    "InputError",
    "InvalidSolventError",
    "NotCompletedError",
    "NotStartedError",
    "SettingConflictError",
]
