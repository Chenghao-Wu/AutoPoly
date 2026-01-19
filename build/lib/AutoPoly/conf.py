#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration Module for AutoPoly Package

This module contains configuration settings for the AutoPoly package including
logging configuration and output path settings.

The configuration includes:
- Output path settings for generated files
- Logging level configurations
- File output settings
- Resource limits to prevent DoS attacks

Configuration can be modified by changing the values in the LOG dictionary
or by setting the OUT_PATH variable.
"""
import logging
import os

# Output path configuration
OUT_PATH = os.path.join(os.path.expanduser("~"))

# Logging configuration dictionary
LOG = {
    'ROOT_LEVEL': logging.INFO,      # Root logger level
    'CONSOLE_LEVEL': logging.INFO,   # Console output level
    'FILE_LEVEL': logging.INFO,      # File output level
    'TO_FILE': False                 # Whether to log to file (True/False)
}

# Resource limits to prevent DoS attacks
# These limits prevent resource exhaustion from excessively large values
MAX_DOP = 10000  # Maximum degree of polymerization
MAX_SEQUENCE_LENGTH = 10000  # Maximum sequence length
MAX_UNIQUE_MONOMERS = 100  # Maximum unique monomer types