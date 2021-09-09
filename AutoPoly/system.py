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

class System(object):
    def __init__(self,out=None):
        self.Folder    =   out
        self.MadeFolder     =   False
        
        self.create_FolderPath()
    @property
    def get_FolderPath(self):
        Folder =   self.Folder
        SystemPath   =   os.getcwd()+"/"+Folder
        return SystemPath

    def create_FolderPath(self):
        self.MadeFolder=True

        FolderPath=self.get_FolderPath
        path=Path(FolderPath)
        if path.exists():
            response = input(FolderPath+" folder exist, delete and make new?(y/n) ")
            if response[0] == "y":
                proc = subprocess.Popen(['/bin/bash'], shell=True,stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                stdout = proc.communicate(("rm -r "+ FolderPath).encode())
                logger.info(' '.join(["removing "+FolderPath]))
                time.sleep(3)
                path.mkdir(parents=True, exist_ok=True)

            elif response[0] == "n":
                logger.info(' '.join(["EXIT : "+FolderPath+" has already existed"]))
                sys.exit()
        else:
            path.mkdir(parents=True, exist_ok=True)

    