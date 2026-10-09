import logging
from collections.abc import Callable

import numpy as np

try:
    from gfnff import GFNFFCalculator
except ImportError:
    GFNFFCalculator = None

from stko._internal.internal_types import MoleculeT
from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.utilities.exceptions import WrapperNotInstalledError

logger = logging.getLogger(__name__)

# The gfnff library works in atomic units.
_BOHR = 0.52917721067  # Å
_HARTREE = 27.211386245988  # eV

# Square-planar restraints: 4-coordinate Pd and Pt keep each cis angle within
# 20 degrees of 90 and each trans angle within 20 degrees of 180.
_SQUARE_PLANAR_METALS = {46, 78}
_SQUARE_PLANAR_TOLERANCE = 20.0  # degrees
_SQUARE_PLANAR_K = 10.0  # eV/rad²


class GFNFFTopo(Optimizer):
    """Uses GFN-FF with the molecule's own bonds as its topology.

    In xTB, as in :class:`.XTBFF`, GFN-FF sets up its topology from the
    coordinates and stores it in the binary ``gfnff_topo`` file, which
    cannot be edited. :class:`GFNFFTopo` instead builds this topology from
    the bonds of the :class:`stk.Molecule` (the bond graph). Stretched
    bonds, such as those between building blocks of a
    :class:`stk.ConstructedMolecule` or metal-ligand bonds, are therefore
    pulled in rather than missed or broken, and atoms that start too close
    cannot form wrong bonds (e.g. hydrogens moving to the wrong atom).

    GFN-FF energies and forces are computed by the standalone ``gfnff``
    library, without an xTB executable. As the library has no optimiser,
    the structure is minimised with L-BFGS.

    See Also:
        * GFN-FF: S. Spicher, S. Grimme, *Angew. Chem. Int. Ed.* **2020**,
          59, 15665, https://doi.org/10.1002/anie.202004239
        * gfnff library: https://github.com/pprcht/gfnff
        * L-BFGS: D. C. Liu, J. Nocedal, *Math. Program.* **1989**, 45,
          503, https://doi.org/10.1007/BF01589116

    Parameters:
        version:
            The GFN-FF parametrisation, e.g. ``"conformer2020"`` (full
            GFN-FF, in which bonds cannot dissociate) or ``"harmonic2020"``
            (a simplified, harmonic GFN-FF built from the bond graph, which
            pulls long bonds in directly). ``None`` uses the library default.

        use_bond_graph:
            If ``True``, the bonds of the molecule are used as the GFN-FF
            topology. If ``False``, GFN-FF works out the bonds from the
            coordinates, as in :class:`.XTBFF`.

        charge:
            Formal molecular charge.

        solvent:
            Implicit solvent name, e.g. ``"h2o"``. An empty string runs in
            vacuum.

        fmax:
            Convergence criterion: the largest force on any atom, in eV/Å.

        max_steps:
            The maximum number of optimisation steps.

        max_step:
            The largest distance any atom may move in one step, in Å.

        square_planar:
            If ``True``, square-planar Pd and Pt are kept square planar with
            flat-bottom angle restraints (cis within 20° of 90°, trans
            within 20° of 180°). Without this, they relax towards
            tetrahedral, as GFN-FF, a force field, does not explicitly
            calculate their electronic structure. Off by default.

    Notes:
        Requires ``gfnff`` 0.3.0 or later (``pip install gfnff``). Wheels
        exist for Linux and Apple-silicon macOS (native arm64 Python).
        Elsewhere, pip builds it from source, which needs a Fortran
        compiler.

        While no pre-optimiser is required, a reasonable starting structure
        (e.g. from :class:`stk.MCHammer`) still helps, as the optimisation
        only finds the nearest local minimum. A poor start, such as
        overlapping building blocks, can end in a strained or tangled
        structure.

        For large molecules, ``gfnff`` can crash if the thread stack is too
        small. If so, set the ``OMP_STACKSIZE`` environment variable (e.g.
        to ``"4000M"``) before ``gfnff`` is imported.

        This optimiser was implemented in stko by Paula C. P. Teeuwen, using
        the ``gfnff`` library. It is adapted from MOCCA, where it was first
        implemented for modelling supramolecular cages (P. C. P. Teeuwen,
        MOCCA, version 1.0, Zenodo, 2026,
        https://doi.org/10.5281/zenodo.21144486).

    Examples:
        For a :class:`stk.ConstructedMolecule`, the bonds made during
        construction are often much too long. A harmonic GFN-FF stage pulls
        these bonds in, after which the full force field finishes the
        optimisation:

        .. code-block:: python

            import stk
            import stko

            bb1 = stk.BuildingBlock('NCCNCCN', [stk.PrimaryAminoFactory()])
            bb2 = stk.BuildingBlock('O=CCCC=O', [stk.AldehydeFactory()])
            polymer = stk.ConstructedMolecule(
                stk.polymer.Linear(
                    building_blocks=(bb1, bb2),
                    repeating_unit="AB",
                    orientations=[0, 0],
                    num_repeating_units=1
                )
            )

            gfnff = stko.OptimizerSequence(
                stko.GFNFFTopo(version='harmonic2020'),
                stko.GFNFFTopo(version='conformer2020'),
            )
            polymer = gfnff.optimize(polymer)

        For a charged molecule, such as a metal-organic cage, set the total
        charge. This example optimises a Pd2L4 lantern (4+) in implicit
        DMSO, with `square_planar` keeping the Pd(II) centres square planar:

        .. code-block:: python

            import numpy as np

            palladium = stk.BuildingBlock(
                smiles='[Pd+2]',
                functional_groups=(
                    stk.SingleAtom(stk.Pd(0, charge=2)) for _ in range(4)
                ),
                position_matrix=np.array([[0.0, 0.0, 0.0]]),
            )
            ligand = stk.BuildingBlock(
                smiles='C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1',
                functional_groups=[
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[#6]~[#7X2]~[#6]',
                        bonders=(1, ),
                        deleters=(),
                    ),
                ],
            )
            lantern = stk.ConstructedMolecule(
                stk.cage.M2L4Lantern(building_blocks=(palladium, ligand)),
            )

            gfnff = stko.OptimizerSequence(
                stko.GFNFFTopo(
                    version='harmonic2020',
                    charge=4,
                    solvent='dmso',
                    square_planar=True,
                ),
                stko.GFNFFTopo(
                    version='conformer2020',
                    charge=4,
                    solvent='dmso',
                    square_planar=True,
                ),
            )
            lantern = gfnff.optimize(lantern)

    """

    def __init__(  # noqa: PLR0913
        self,
        *,
        version: str | None = "conformer2020",
        use_bond_graph: bool = True,
        charge: int = 0,
        solvent: str = "",
        fmax: float = 0.05,
        max_steps: int = 5000,
        max_step: float = 0.2,
        square_planar: bool = False,
    ) -> None:
        if GFNFFCalculator is None:
            msg = (
                "GFNFFTopo needs the gfnff library, version 0.3.0 or later "
                "(pip install gfnff)."
            )
            raise WrapperNotInstalledError(msg)

        self._version = version
        self._use_bond_graph = use_bond_graph
        self._charge = charge
        self._solvent = solvent
        self._fmax = fmax
        self._max_steps = max_steps
        self._max_step = max_step
        self._square_planar = square_planar

    def _bond_matrix(self, mol: MoleculeT) -> np.ndarray:
        """Every bond of `mol`, including order-0 (e.g. dative) bonds."""
        num_atoms = mol.get_num_atoms()
        bond_matrix = np.zeros((num_atoms, num_atoms), dtype=np.int32)
        for bond in mol.get_bonds():
            i = bond.get_atom1().get_id()
            j = bond.get_atom2().get_id()
            bond_matrix[i, j] = bond_matrix[j, i] = 1
        return bond_matrix

    def _square_planar_restraints(
        self,
        mol: MoleculeT,
    ) -> list[tuple[int, int, int, float, float, float]]:
        """Angle restraints keeping 4-coordinate Pd and Pt square planar."""
        neighbours: dict[int, list[int]] = {}
        for bond in mol.get_bonds():
            i = bond.get_atom1().get_id()
            j = bond.get_atom2().get_id()
            neighbours.setdefault(i, []).append(j)
            neighbours.setdefault(j, []).append(i)
        numbers = [atom.get_atomic_number() for atom in mol.get_atoms()]
        positions = mol.get_position_matrix()
        restraints = []
        for centre, ligands in neighbours.items():
            if numbers[centre] not in _SQUARE_PLANAR_METALS:
                continue
            if len(ligands) != 4:  # noqa: PLR2004
                continue
            pairs = []
            for a, atom1 in enumerate(ligands):
                for atom3 in ligands[a + 1 :]:
                    u = positions[atom1] - positions[centre]
                    v = positions[atom3] - positions[centre]
                    cosine = u @ v / (np.linalg.norm(u) * np.linalg.norm(v))
                    angle = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
                    pairs.append((atom1, atom3, angle > 135))  # noqa: PLR2004
            if sum(trans for _, _, trans in pairs) != 2:  # noqa: PLR2004
                msg = (
                    f"Atom {centre} does not start square planar, so it is "
                    "not restrained."
                )
                logger.warning(msg)
                continue
            restraints += [
                (
                    atom1,
                    centre,
                    atom3,
                    180.0 if trans else 90.0,
                    _SQUARE_PLANAR_TOLERANCE,
                    _SQUARE_PLANAR_K,
                )
                for atom1, atom3, trans in pairs
            ]
        return restraints

    def _restraint_forces(
        self,
        positions: np.ndarray,
        angle_restraints: list[tuple[int, int, int, float, float, float]],
    ) -> np.ndarray:
        """Forces (eV/Å) of the angle restraints, positions in Å."""
        forces = np.zeros_like(positions)
        for atom1, centre, atom3, angle, tolerance, k in angle_restraints:
            u = positions[atom1] - positions[centre]
            v = positions[atom3] - positions[centre]
            ru, rv = np.linalg.norm(u), np.linalg.norm(v)
            cosine = np.clip(u @ v / (ru * rv), -1.0, 1.0)
            difference = np.arccos(cosine) - np.radians(angle)
            excess = np.sign(difference) * max(
                abs(difference) - np.radians(tolerance), 0.0
            )
            if excess == 0:
                continue
            # dE/dcos = dE/dtheta * dtheta/dcos. The window keeps active
            # angles away from 0 and 180 degrees.
            sine = np.sqrt(max(1.0 - cosine**2, 1e-12))
            de_dcos = -k * excess / sine
            gradient1 = de_dcos * (v / (ru * rv) - cosine * u / ru**2)
            gradient3 = de_dcos * (u / (ru * rv) - cosine * v / rv**2)
            forces[atom1] -= gradient1
            forces[atom3] -= gradient3
            forces[centre] += gradient1 + gradient3
        return forces

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimise `mol`.

        Parameters:
            mol:
                The molecule to be optimised.

        Returns:
            The optimised molecule.

        """
        numbers = np.array(
            [atom.get_atomic_number() for atom in mol.get_atoms()],
            dtype=np.int32,
        )
        positions = np.ascontiguousarray(mol.get_position_matrix() / _BOHR)
        calculator = GFNFFCalculator(
            numbers,
            positions,
            charge=self._charge,
            solvent=self._solvent,
            version=self._version,
            bond_matrix=(
                self._bond_matrix(mol) if self._use_bond_graph else None
            ),
        )

        angle_restraints = (
            self._square_planar_restraints(mol) if self._square_planar else []
        )

        def forces(positions: np.ndarray) -> np.ndarray:
            """GFN-FF and restraint forces in eV/Å, positions in Å."""
            positions_au = np.ascontiguousarray(positions / _BOHR)
            _, gradient, _ = calculator.singlepoint(numbers, positions_au)
            return -gradient * _HARTREE / _BOHR + self._restraint_forces(
                positions, angle_restraints
            )

        try:
            positions, steps, largest_force = lbfgs_minimize(
                forces,
                mol.get_position_matrix(),
                fmax=self._fmax,
                max_steps=self._max_steps,
                max_step=self._max_step,
            )
        finally:
            calculator.deallocate()

        if largest_force > self._fmax:
            msg = (
                f"GFNFFTopo did not converge to fmax={self._fmax} eV/Å "
                f"within {steps} steps (largest force "
                f"{largest_force:.3f} eV/Å)."
            )
            logger.warning(msg)

        return mol.with_position_matrix(positions)


def lbfgs_minimize(  # noqa: PLR0913
    forces: Callable[[np.ndarray], np.ndarray],
    positions: np.ndarray,
    *,
    fmax: float,
    max_steps: int,
    max_step: float,
    memory: int = 100,
    alpha: float = 70.0,
) -> tuple[np.ndarray, int, float]:
    """Minimise with L-BFGS, without a line search.

    The search direction comes from the L-BFGS two-loop recursion (D. C.
    Liu, J. Nocedal, *Math. Program.* **1989**, 45, 503,
    https://doi.org/10.1007/BF01589116) over the last `memory` steps,
    starting from a diagonal Hessian `alpha` (eV/Å²). Instead of a line
    search, each step is scaled down so that no atom moves more than
    `max_step` (Å), as a line search copes poorly with GFN-FF's not
    entirely smooth energy. Updates with non-positive curvature are
    skipped.

    The approach without a line search and the defaults follow ASE's
    ``LBFGS`` (A. H. Larsen et al., *J. Phys.: Condens. Matter* **2017**,
    29, 273002, https://doi.org/10.1088/1361-648X/aa680e). The code is
    written independently.

    Parameters:
        forces:
            Forces (eV/Å) for positions (Å), both of shape ``(n, 3)``.

        positions:
            Starting positions (Å).

        fmax:
            Converged once the largest force on any atom is below this.

        max_steps:
            The maximum number of steps.

        max_step:
            The largest distance any atom may move in one step (Å).

        memory:
            The number of previous steps used.

        alpha:
            Initial guess of the Hessian's diagonal (eV/Å²).

    Returns:
        The final positions, the number of steps taken and the largest
        force on any atom.

    """
    s_list: list[np.ndarray] = []
    y_list: list[np.ndarray] = []
    rho_list: list[float] = []
    f = forces(positions)
    steps = 0
    while True:
        largest_force = float(np.linalg.norm(f, axis=1).max())
        if largest_force < fmax or steps >= max_steps:
            return positions, steps, largest_force

        q = -f.ravel()
        a = []
        for s, y, rho in zip(
            reversed(s_list),
            reversed(y_list),
            reversed(rho_list),
            strict=True,
        ):
            a.append(rho * s @ q)
            q = q - a[-1] * y
        z = q / alpha
        for s, y, rho, a_i in zip(
            s_list, y_list, rho_list, reversed(a), strict=True
        ):
            z = z + s * (a_i - rho * y @ z)
        step = -z.reshape(-1, 3)

        longest = np.linalg.norm(step, axis=1).max()
        if longest > max_step:
            step *= max_step / longest

        new_positions = positions + step
        new_f = forces(new_positions)
        s = (new_positions - positions).ravel()
        y = (f - new_f).ravel()
        if y @ s > 0:  # keep the curvature information positive
            s_list.append(s)
            y_list.append(y)
            rho_list.append(1.0 / (y @ s))
            if len(s_list) > memory:
                del s_list[0], y_list[0], rho_list[0]
        positions, f = new_positions, new_f
        steps += 1
