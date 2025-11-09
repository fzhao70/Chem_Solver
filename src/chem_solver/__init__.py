"""
SMVGEAR-like Chemical Solver Package

A chemical kinetics solver for atmospheric chemistry simulations.
"""

from .solver import ChemicalSolver
from .mechanism import ChemicalMechanism

__version__ = "0.1.0"
__all__ = ["ChemicalSolver", "ChemicalMechanism"]
