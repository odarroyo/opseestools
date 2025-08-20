"""
opseestools - A Python library to streamline OpenSeesPy workflows

This package provides tools for structural analysis using OpenSeesPy.

Modules:
    analisis: 2D analysis functions 
    analisis3D: 3D analysis functions
    utilidades: Utility functions for model building
    Lib_frag: Fragility and vulnerability analysis functions
"""

__version__ = "0.7.0"
__author__ = "opseestools contributors"

# Try to import main modules for easier access
# Handle missing dependencies gracefully for documentation building
try:
    from . import analisis
    from . import analisis3D  
    from . import utilidades
    from . import Lib_frag
except ImportError:
    # This allows Sphinx autodoc to work even when dependencies are missing
    pass