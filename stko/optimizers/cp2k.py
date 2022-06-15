"""This module defines an stko interface to CP2K

http://www.cp2k.org

Author: Steven Bennett <s.bennett18@imperial.ac.uk>
"""
import logging
import os
import re
import shutil
import subprocess as sp
from uuid import uuid4

import numpy as np

from rdkit.Chem import AllChem as rdkit

from .optimizers import Optimizer

logger = logging.getLogger(__name__)


class CP2K(Optimizer):
    """
    stko optimizer for CP2K.

    CP2K is a program to perform atomistic and molecular simulations of
    solid state, liquid, molecular, and biological systems.
    It provides a general
    framework for different methods such as e.g., density functional
    theory (DFT) using a mixed Gaussian and plane waves approach (GPW)
    and classical pair and many-body potentials.
    """

    def __init__(
        self,
        cp2k_path,
        input_template,
        output_dir=None,
        job_name=None,
        run_type="GEO_OPT",
    ):
        """
        Initialize CP2K optimizer.

        Parameters
        ----------
        cp2k_path : :class:`str`
            Path to CP2K executable file.

        input_template : :class:`str`
            Path to a CP2K input template to read.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        job_name : :class:`str`, optional
            The name of the job, if ``None`` then :func:`uuid.uuid4` is
            used.

        run_type : :class:`str`, optional
            The type of calculation to perform. `GEO_OPT` is the
            default.

        """
        self._cp2k_path = cp2k_path
        self._input_template = input_template
        self._output_dir = output_dir
        self._job_name = job_name
        self._run_type = run_type

    def _write_input_file(self, mol, in_file):
        """
        Write CP2K input file.

        Parameters
        ----------
        mol : :class:`stk.Molecule`
            Molecule to optimise.

        in_file : :class:`str`
            Path to write the CP2K input file.
        """
        with open(self._input_template, "r") as f:
            input_template = f.read()

        # Replace the global section of the input file.
        global_section = self._get_global_section()
        input_template = re.sub(
            r"&GLOBAL.*?&END GLOBAL",
            global_section,
            input_template,
            flags=re.DOTALL,
        )

        # Replace the coordinate section of input file
        coordinate_section = self._get_coord_section(mol)
        # Replaces between the first occurrence of "&COORD" and
        # "&END COORD"
        input_template = re.sub(
            r"&COORD.*&END COORD",
            coordinate_section,
            input_template,
            flags=re.DOTALL,
        )
        # Write the input file
        with open(in_file, "w") as f:
            f.write(input_template)
        # Replace the cell section of the input file
        cell_section = self._get_cell_section(mol)
        # Replaces between the first occurrence of "&CELL"
        # and "&END CELL"
        input_template = re.sub(
            r"&CELL.*&END CELL",
            cell_section,
            input_template,
            flags=re.DOTALL,
        )
        # Write the modified input string
        with open(in_file, "w") as f:
            f.write(input_template)

    def _get_global_section(self):
        """
        Return the global section of the CP2K input file.

        Returns
        -------
        :class:`str`
            Title section of the CP2K input file.
        """
        if self._job_name is None:
            self._job_name = str(uuid4().int)
        title_section = "&GLOBAL\n"
        title_section += f"  PROJECT_NAME {self._job_name}\n"
        title_section += f"  RUN_TYPE {self._run_type}\n"
        title_section += "&END GLOBAL"
        return title_section

    def _get_cell_section(self, mol):
        """
        Return the cell section of the CP2K input file.

        Parameters
        ----------
        mol : :class:`stk.Molecule`
            Molecule to write.

        Returns
        -------
        :class:`str`
            Cell section of the CP2K input file.
        """
        cell_section = "&CELL\n"
        # Find the furthest cartesian coordinate from the origin in all
        # directions.
        position_matrix = mol.get_position_matrix()
        # For a single molecule, ensure the box size is equal to the
        # furthest coordinate from the origin in doubled.
        furthest_coord = round(np.max(np.max(position_matrix, axis=0)), 1)
        cell_pos = int(furthest_coord + 8)
        cell_section += f"ABC {cell_pos} " f"{cell_pos} {cell_pos}\n"
        cell_section += "PERIODIC NONE\n"
        cell_section += "&END CELL"
        return cell_section

    def _get_coord_section(self, mol):
        """
        Return the coordinate section of the CP2K input file.

        Parameters
        ----------
        mol : :class:`stko.molecule.Molecule`
            Molecule to write.

        Returns
        -------
        :class:`str`
            Coordinate section of the CP2K input file.
        """
        coord_section = "&COORD\n"
        for atom in mol.get_atoms():
            atom_symbol = rdkit.Atom(atom.get_atomic_number()).GetSymbol()
            position = mol.get_centroid(atom_ids=atom.get_id())
            coord_section += (
                f"{atom_symbol}    {round(position[0], 5)}    "
                f"{round(position[1], 5)}"
                f"      {round(position[2], 5)}\n"
            )
        coord_section += "&END COORD"
        return coord_section

    def optimize(self, mol, run_name=str(uuid4().int)):
        """
        Optimize a molecule.

        Parameters
        ----------
        mol : :class:`stko.molecule.Molecule`
            Molecule to optimise.

        run_name : :class:`str`, optional
            The name of the run files.

        Returns
        -------
        :class:`stko.molecule.Molecule`
            Optimized molecule.
        """
        # TODO Write input file
        if self._output_dir is None:
            output_dir = str(uuid4().int)
        else:
            output_dir = self._output_dir
        if os.path.exists(output_dir):
            # Add overwrite warning to documentation
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        in_file = f"{run_name}.inp"
        log_file = f"{run_name}.log"

        if self._input_template is not None:
            input_template = self._input_template
        else:
            raise NotImplementedError(
                "No input template specified."
                "Support for generating input files from scratch"
                " is not yet implemented."
            )
        try:
            self._write_input_file(
                input_template=input_template, mol=mol, in_file=in_file
            )
            self._run_cp2k(in_file, log_file)
            # Update from output file.
            # TODO: Check output xyz file name
            mol = mol.with_structure_from_file(
                self._extract_positions()  # Add input template name
            )

        finally:
            os.chdir(init_dir)

        return mol

    def _extract_positions(self, output_xyz):
        """
        Extract the positions from the CP2K output file.

        Parameters
        ----------
        output_xyz : :class:`str`
            Path to the CP2K output xyz file.

        Returns
        -------
        :class:`numpy.ndarray`
            Positions of the atoms.
        """
        with open(output_xyz, "r") as f:
            lines = f.readlines()
        # Find index of final line with only numbers and spaces
        final_xyz_line = [
            lines.index(line, i)
            for i, line in enumerate(lines)
            if re.match(r"^[0-9\s]+$", line)
        ][-1]
        xyz_string = [lines[i] for i in range(final_xyz_line, len(lines))]
        final_xyz = output_xyz.replace(".xyz", "_final.xyz")
        with open(final_xyz, "w") as f:
            f.write("".join(xyz_string))
        return final_xyz

    def _run_cp2k(self, in_file, out_file):
        cmd = f"{self._cp2k_path} -i {in_file} -o {out_file}"
        with open(out_file, "w") as f:
            # Call will hold the the program until execution completion
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell required to run complex arguments.
                shell=True,
            )
