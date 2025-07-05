#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Logging Module for AutoPoly Package

This module provides logging functionality for the AutoPoly package.
It sets up a centralized logging system with configurable output levels
for both console and file logging.

The logging system includes:
- Configurable log levels (DEBUG, INFO, WARNING, ERROR, CRITICAL)
- Console output with formatted timestamps
- Optional file output for persistent logging
- Centralized configuration through the conf module

Usage:
    from .logger import setup_logger
    logger = setup_logger()
    logger.info("Your message here")
"""
import os
import logging
from .conf import LOG, OUT_PATH


def setup_logger(to_file: bool = None) -> logging.Logger:
    """
    Set up and configure the logger for the AutoPoly package.
    
    This function creates a logger instance with the following features:
    - Configurable log levels for root, console, and file handlers
    - Formatted output with timestamps
    - Optional file logging to analysis_log.log
    - Console output for immediate feedback
    
    Args:
        to_file (bool, optional): Whether to enable file logging.
                                 If None, uses the default from conf.LOG.
                                 Defaults to None.
    
    Returns:
        logging.Logger: Configured logger instance
    
    Example:
        >>> logger = setup_logger()
        >>> logger.info("Starting polymer generation...")
        2024-01-15 10:30:45 - INFO - Starting polymer generation...
    """
    # Use default setting if to_file is not specified
    if to_file is None:
        to_file = LOG['TO_FILE']
    
    logger = logging.getLogger(__name__)
    logger.setLevel(LOG['ROOT_LEVEL'])
    
    # Create formatter with timestamp
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        '%Y-%m-%d %H:%M:%S'
    )

    # Add file handler if requested
    if to_file:
        file_handler = logging.FileHandler(
            os.path.join(OUT_PATH, 'analysis_log.log'), 
            mode='w'
        )
        file_handler.setLevel(LOG['FILE_LEVEL'])
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # Add console handler
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(LOG['CONSOLE_LEVEL'])
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    return logger