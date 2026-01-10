#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File management utilities for the polymerization package.

This module provides file and directory management functions for the
polymerization workflow, including directory creation, file movement,
and settings file modification.
"""

import os
from pathlib import Path
import shutil
from .system import logger
from .exceptions import WorkflowError


def create_working_directory(system, name):
    """
    Create and manage the working directory structure for the polymerization.

    This function sets up the directory structure needed for the polymerization
    process. It creates the main project directory and the moltemplate
    subdirectory where all intermediate files will be stored.

    The function handles existing directories by prompting the user to either
    delete and recreate them or choose a different project name.

    Args:
        system: System object with get_folder_path() method
        name: Name of the polymerization project

    Returns:
        Path object to the moltemplate directory

    Raises:
        SystemExit: If user chooses not to overwrite existing directory
    """
    base_path = Path(system.get_folder_path())
    polymer_path = base_path / name
    moltemplate_path = polymer_path / "moltemplate"

    # Check if base directory exists
    if polymer_path.exists():
        response = input(f"{polymer_path} folder exists, delete and make new?(y/n) ")
        if response.lower() == 'y':
            logger.info(f"removing {polymer_path}")
            shutil.rmtree(polymer_path)
        else:
            raise WorkflowError(
                "Directory exists. Please remove the existing folder or choose a different name."
            )

    # Create directory structure
    moltemplate_path.mkdir(parents=True, exist_ok=True)
    return moltemplate_path


def get_rid_of_lj_cut_coul_long(path_cwd):
    """
    Removes lj/cut/coul/long from the settings file.

    This function reads the system.in.settings file, removes any
    lj/cut/coul/long pair_coeff entries, and rewrites the file.

    Args:
        path_cwd: Current working directory path

    Raises:
        SystemExit: If system.in.settings file does not exist
    """
    in_ = Path(path_cwd) / "system.in.settings"
    out = Path(path_cwd) / "tmp.data"

    in_path = Path(in_)
    if not in_path.is_file():
        raise WorkflowError(f"system.in.settings does not exist: {in_}")

    # Use context manager to ensure file is properly closed even if exception occurs
    with open(in_, 'r') as read_f, open(out, "w") as write_f:
        while True:
            line = read_f.readline()
            if line.strip() == "":
                write_f.write("\n")
            elif line.strip().split()[0] == "pair_coeff":
                #write_f.write("    pair_coeff ")
                space_i = 0
                for ii in line.split():
                    if ii == "lj/cut/coul/long":
                        continue
                    else:
                        if space_i == 0:
                            write_f.write("    ")
                            space_i = space_i + 1
                        write_f.write(ii + " ")
                write_f.write("\n")
            else:
                write_f.write(line)

            if not line:
                break
    # Use atomic Path.replace() instead of unsafe shell command
    Path(out).replace(in_)


def mv_files(path_cwd):
    """
    Moves generated files to the appropriate directories.

    This function organizes the output files from moltemplate by:
    1. Copying system files to the parent directory
    2. Moving system input/output files to the output directory
    3. Moving .lt and .prm files to the input directory

    Args:
        path_cwd: Current working directory path (moltemplate directory)

    Raises:
        SystemExit: If an error occurs while moving files
    """
    try:
        # Define paths
        moltemplate_dir = Path(path_cwd)
        parent_dir = moltemplate_dir.parent

        # Create output and input directories if they don't exist
        output_dir = parent_dir / "output"
        input_dir = parent_dir / "input"
        output_dir.mkdir(exist_ok=True)
        input_dir.mkdir(exist_ok=True)

        # Copy data files to parent directory
        for file in ["system.data", "system.in.charges", "system.in.settings", "system.in", "system.in.init"]:
            if (moltemplate_dir / file).exists():
                shutil.copy2(moltemplate_dir / file, parent_dir)

        # Move files to output directory
        for pattern in ["system.in*", "system*data", "output_ttree"]:
            for file in moltemplate_dir.glob(pattern):
                shutil.move(str(file), str(output_dir))

        # Move files to input directory
        for pattern in ["*.lt", "*.prm"]:
            for file in moltemplate_dir.glob(pattern):
                shutil.move(str(file), str(input_dir))

    except Exception as e:
        raise WorkflowError(f"Error moving files: {str(e)}") from e
