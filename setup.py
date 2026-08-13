# Package metadata lives in pyproject.toml (PEP 621).
# This shim only exists so legacy tooling that requires setup.py keeps working.
from setuptools import setup

setup()
