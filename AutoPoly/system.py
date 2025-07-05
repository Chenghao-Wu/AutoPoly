#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
System Management Module

This module provides the System class for managing file paths, directories,
and system operations within the AutoPoly package.

The System class handles:
- Output directory creation and management
- Path resolution for various file operations
- Integration with the logging system

Created on Fri Dec 21 12:19:08 2018
@author: zwu
"""
import os
import sys
from pathlib import Path 
import subprocess
import time
import shutil

import subprocess

from .logger import setup_logger
logger = setup_logger()

class System:
    """
    System management class for AutoPoly package.
    
    This class handles file system operations, path management, and provides
    utilities for creating and managing output directories for polymer simulations.
    
    Attributes:
        out (str): Output directory name for the current simulation
        get_FolderPath (str): Full path to the output directory
    """
    
    def __init__(self, out: str = None) -> None:
        """
        Initialize the System class.
        
        Args:
            out (str, optional): Output directory name. If None, uses current
                               working directory. Defaults to None.
        """
        self.out = out
        self.get_FolderPath = os.getcwd() + "/" + self.out if self.out else os.getcwd()
        
        # Create output directory if specified
        if self.out:
            self._create_output_directory()
    
    def _create_output_directory(self) -> None:
        """
        Create the output directory if it doesn't exist.
        
        This method creates the output directory specified in the constructor.
        If the directory already exists, it logs a warning but continues.
        """
        try:
            Path(self.get_FolderPath).mkdir(parents=True, exist_ok=True)
            logger.info(f"Output directory created/verified: {self.get_FolderPath}")
        except Exception as e:
            logger.error(f"Failed to create output directory {self.get_FolderPath}: {e}")
            raise
    
    def get_output_path(self) -> str:
        """
        Get the full path to the output directory.
        
        Returns:
            str: Full path to the output directory
        """
        return self.get_FolderPath
    
    def change_output_directory(self, new_out: str) -> None:
        """
        Change the output directory to a new location.
        
        Args:
            new_out (str): New output directory name
        """
        self.out = new_out
        self.get_FolderPath = os.getcwd() + "/" + self.out
        self._create_output_directory()
        logger.info(f"Output directory changed to: {self.get_FolderPath}")
    
    def cleanup_output_directory(self) -> None:
        """
        Remove the output directory and all its contents.
        
        This method should be used with caution as it permanently deletes
        all files in the output directory.
        """
        if os.path.exists(self.get_FolderPath):
            try:
                shutil.rmtree(self.get_FolderPath)
                logger.info(f"Cleaned up output directory: {self.get_FolderPath}")
            except Exception as e:
                logger.error(f"Failed to cleanup output directory: {e}")
                raise
