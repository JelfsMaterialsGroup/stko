import logging
import pathlib
from typing import Protocol

import stk

from stko._internal.internal_types import MoleculeT

logger = logging.getLogger(__name__)


class Optimizer(Protocol):
    """A base class for optimizers.

    Optimizers are objects used to optimize molecules. Each optimizer is
    initialized with some settings and can optimize a molecule
    with :meth:`~.Optimizer.optimize`.

    .. testcode:: optimiser-example

        import stk
        import stko

        mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
        mmff = stko.MMFF()
        mol = mmff.optimize(mol)

        # Optimizers also work with ConstructedMolecule objects.
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(mol, ),
                repeating_unit='A',
                num_repeating_units=3,
            )
        )
        etkdg = stko.ETKDG()
        polymer = etkdg.optimize(polymer)

    Sometimes it is desirable to chain multiple optimizations, one after
    another. For example, before running an optimization, it may be
    desirable to embed a molecule first, to generate an initial structure.
    :class:`.OptimizerSequence` may be used for this.

    .. testcode:: optimiser-example

        # Create a new optimizer which chains the previously defined
        # mmff and etkdg optimizers.
        optimizer_sequence = stko.OptimizerSequence(etkdg, mmff)

        # Run each optimizer in sequence.
        polymer = optimizer_sequence.optimize(polymer)

    .. _`adding optimizers`:

    Making New Optimizers:
    ----------------------

    New optimizers can be made by simply making a class which defines a
    :meth:`~.Optimizer.optimize` method. The method must take 1
    mandatory `mol` parameter. :meth:`~.Optimizer.optimize` will take the
    `mol` and change its structure in whatever way it likes. Beyond this
    there are no requirements. New optimizers can be added into the
    :mod:`.optimizers` submodule or into a new submodule.

    """

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol: The molecule to be optimized.

        Returns:
            The optimized molecule.

        """
        raise NotImplementedError


class NullOptimizer(Optimizer):
    """Applies no optimizer."""

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        return mol


class OptimizerSequence(Optimizer):
    """Applies optimizers in sequence.

    Parameters:
        optimizers:
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

    Examples:
        Let's say we want to embed a molecule with ETKDG first and then
        minimize it with the MMFF force field.

        .. testcode:: optimiser-sequence

            import stk
            import stko

            mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
            optimizer = stko.OptimizerSequence(stko.ETKDG(), stko.MMFF())
            mol = optimizer.optimize(mol)

    """

    def __init__(self, *optimizers: Optimizer) -> None:
        self._optimizers = optimizers

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        for optimizer in self._optimizers:
            cls_name = optimizer.__class__.__name__
            msg = f'Using {cls_name} on "{mol}".'
            logger.info(msg)
            mol = optimizer.optimize(mol)

        return mol


class OptWriterSequence(Optimizer):
    """Applies optimizers in sequence and writes each step to avoid reruns.

    Parameters:
        optimizers:
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

    Examples:
        Let's say we want to embed a molecule with ETKDG first and then
        minimize it with the MMFF force field.

        .. testcode:: optimiser-writer-sequence

            import stk
            import stko
            import pathlib

            output_directory = pathlib.Path('output_path')
            output_directory.mkdir(exist_ok=True)

            mol = stk.BuildingBlock('NCCCN', [stk.PrimaryAminoFactory()])
            optimizer = stko.OptWriterSequence(
                optimizers={
                    "etkdg": stko.ETKDG(),
                    "mmff": stko.MMFF(),
                },
                writer=stk.MolWriter(),
                output_directory=output_directory,
            )
            mol = optimizer.optimize(mol)

        .. testcleanup:: optimiser-writer-sequence

            import shutil

            shutil.rmtree('output_path')

    """

    def __init__(
        self,
        optimizers: dict[str, Optimizer],
        writer: stk.XyzWriter
        | stk.MolWriter
        | stk.PdbWriter
        | stk.TurbomoleWriter,
        output_directory: pathlib.Path,
    ) -> None:
        self._optimizers = optimizers
        self._writer = writer
        self._output_directory = output_directory

    def _suffix_from_writer(self) -> str:
        if isinstance(self._writer, stk.XyzWriter):
            return "xyz"
        if isinstance(self._writer, stk.TurbomoleWriter):
            return "coord"
        if isinstance(self._writer, stk.PdbWriter):
            return "pdb"
        if isinstance(self._writer, stk.MolWriter):
            return "mol"
        return None

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        filesuffix = self._suffix_from_writer()
        for optimizer_name in self._optimizers:
            output_file_name = (
                self._output_directory / f"{optimizer_name}_out.{filesuffix}"
            )
            cls_name = self._optimizers[optimizer_name].__class__.__name__
            if not output_file_name.exists():
                msg = f"Running {cls_name} optimisation."
                logger.info(msg)
                mol = self._optimizers[optimizer_name].optimize(mol)
                self._writer.write(molecule=mol, path=output_file_name)
            else:
                msg = f"Loading {cls_name} optimisation."
                logger.info(msg)
                mol = mol.with_structure_from_file(output_file_name)

        return mol


class TryCatchOptimizer(Optimizer):
    """Try to optimize with a Optimizer, use another on failure.

    Parameters:
        try_optimizer
            The optimizer which is used initially to try and optimize a
            :class:`stk.Molecule`.

        catch_optimizer:
            If `try_optimizer` raises an error, this optimizer is
            run on the :class:`stk.Molecule` instead.


    Examples:
        .. testcode:: try-catch

            import stk
            import stko

            # Create some molecules to optimize.
            mol1 = stk.BuildingBlock('NCCN')
            mol2 = stk.BuildingBlock('CCCCC')
            mol3 = stk.BuildingBlock('O=CCCN')

            # Create an optimizer which may fail.
            uff = stko.UFF()

            # Create a backup optimizer.
            mmff = stko.MMFF()

            # Make an optimizer which tries to run raiser and if that
            # raises an error, will run mmff on the molecule instead.
            try_catch = stko.TryCatchOptimizer(
                try_optimizer=uff,
                catch_optimizer=mmff,
            )

            # Optimize the molecules. In each case if the optimization with
            # UFF fails, MMFF is used to optimize the molecule instead.
            mol1 = try_catch.optimize(mol1)
            mol2 = try_catch.optimize(mol2)
            mol3 = try_catch.optimize(mol3)

    """

    def __init__(
        self,
        try_optimizer: Optimizer,
        catch_optimizer: Optimizer,
    ) -> None:
        self._try_optimizer = try_optimizer
        self._catch_optimizer = catch_optimizer

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        """Optimize `mol`.

        Parameters:
            mol: The molecule to be optimized.

        Returns:
            The molecule to be optimized.

        """
        try:
            return self._try_optimizer.optimize(mol)
        except Exception:
            try_name = self._try_optimizer.__class__.__name__
            catch_name = self._catch_optimizer.__class__.__name__
            msg = f"{try_name} failed, trying {catch_name}."
            logger.exception(msg)
            return self._catch_optimizer.optimize(mol)
