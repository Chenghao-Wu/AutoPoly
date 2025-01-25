#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
    def __init__(self, out: str = None) -> None:
        """Initialize the System class.
        
        Args:
            out (str): Output directory name.
        """
        self.out = out
        self.get_FolderPath = os.getcwd() + "/" + self.out
