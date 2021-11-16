"""This module defines an stko interface to CP2K

http://www.cp2k.org

Author: Steven Bennett <s.bennett18@imperial.ac.uk>
"""
from .optimizers import Optimizer
import subprocess as sp
import logging
import numpy as np
from uuid import uuid4
import shutil
import numpy as np
import os
from rdkit.Chem import AllChem as rdkit
import re

logger = logging.getLogger(__name__)

class CP2K(Optimizer):
    """
    stko optimizer for CP2K.

    CP2K is a program to perform atomistic and molecular simulations of solid
    state, liquid, molecular, and biological systems. It provides a general
    framework for different methods such as e.g., density functional theory
    (DFT) using a mixed Gaussian and plane waves approach (GPW) and classical
    pair and many-body potentials.
    """

    def __init__(
        self,
        cp2k_path,
        input_template,
        output_dir=None,
    ):
        """
        Initialize CP2K optimizer.

        Parameters
        ----------
        cp2k_path : :class:`str`
            Path to CP2K executable file.

        basis_set : :class:`str`
            Basis set to use.

        input_template : :class:`str`
            Path to a CP2K input template to read.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the calculation are written, if ``None`` then
            :func:`uuid.uuid4` is used.
        """
        self._cp2k_path = cp2k_path
        self._input_template = input_template
        self._output_dir = output_dir


    def _write_input_file(self, mol, input_template, in_file):
        """
        Write CP2K input file.

        Parameters
        ----------
        mol : :class:`stk.Molecule`
            Molecule to optimise.

        input_template : :class:`str`
            Path to a CP2K input template to read.

        in_file : :class:`str`
            Path to write the CP2K input file.
        """
        with open(input_template, "r") as f:
            input_template = f.read()
        # Replace the coordinate section of input file
        coordinate_section = self._get_coord_section(mol)
        # Replaces between the first occurrence of "&COORD" and "&END COORD"
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
        # Replaces between the first occurrence of "&CELL" and "&END CELL"
        input_template = re.sub(
            r"&CELL.*&END CELL",
            cell_section,
            input_template,
            flags=re.DOTALL,
        )
        # Write the modified input string
        with open(in_file, "w") as f:
            f.write(input_template)

    def _make_coordinates_positive(self, mol):
        """Translates the centroid to make all atom coordinates positive values"""
        position_matrix = mol.get_position_matrix()
        # Find minimum value of matrix columns
        minimum_pos = np.min(position_matrix, axis=0)
        # Translate centroid to make all coordinates positive
        return mol.with_centroid([abs(coord) for coord in minimum_pos])

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
        cell_section = f"&CELL\n"
        # Find the furthest cartesian coordinate from the origin in all directions
        position_matrix = mol.get_position_matrix()
        # For a single molecule, ensure the box size is large enough so the electron density at the sides is 0.
        furthest_coord = round(np.max(np.max(position_matrix, axis=0)) + 1, 5)
        cell_section += f"ABC {furthest_coord} {furthest_coord} {furthest_coord}\n"
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
        coord_section = f"&COORD\n"
        for atom in mol.get_atoms():
            atom_symbol = rdkit.Atom(atom.get_atomic_number()).GetSymbol()
            position = mol.get_centroid(atom_ids=atom.get_id())
            coord_section += (
                f"{atom_symbol}    {abs(round(position[0], 5))}    "
                f"{abs(round(position[1], 5))}      {abs(round(position[2], 5))}\n"
            )
        coord_section += "&END COORD"
        return coord_section

    def optimize(self, mol):
        """
        Optimize a molecule.

        Parameters
        ----------
        mol : :class:`stko.molecule.Molecule`
            Molecule to optimise.

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

        in_file = "cp2k_input.inp"
        out_file = "cp2k_output.log"

        if self._input_template is not None:
            input_template = self._input_template
        else:
            raise NotImplementedError(
                "No input template specified."
                "Support for generating input files from scratch is not yet implemented."
            )
        try:
            # Make all coordinates positive values
            mol = self._make_coordinates_positive(mol)
            self._write_input_file(
                input_template=input_template,
                mol=mol,
                in_file=in_file
            )
            self._run_cp2k(
                in_file,
                out_file
            )
            # Update from output file.
            mol = mol.with_structure_from_file(self._extract_positions("H2O-pos-1.xyz"))

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
        final_xyz_line = [lines.index(line, i) for i, line in enumerate(lines) if re.match(r"^[0-9\s]+$", line)][-1]
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
