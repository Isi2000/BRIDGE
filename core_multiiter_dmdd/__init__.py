"""
Core Multi-Iterative Dynamic Mode Decomposition with Delay Coordinates

This package implements the Multi-dimensional Higher-Order Dynamic Mode Decomposition (HODMD)
with iterative refinement for complex signal analysis without plotting functionality.

Main classes and functions:
- MultiIterativeHODMD: Main class for multi-iterative HODMD analysis
- dmdd_core: Core DMD with delay coordinates algorithm
- hosvd_function: Higher-order SVD decomposition
- tensor_utils: Utility functions for tensor operations
- reconstruction: Functions for DMD reconstruction
"""

from .multi_iter_hodmd import MultiIterativeHODMD
from .dmdd_core import dmdd_core
from .hosvd_function import hosvd_function
from .tensor_utils import hosvd, tprod, ndim_unfold, ndim_fold, svdtrunc
from .reconstruction import cont_reconst, dmd_reconst, calculate_dmd_mode

__version__ = "1.0.0"
__author__ = "Claude Code"

__all__ = [
    'MultiIterativeHODMD',
    'dmdd_core', 
    'hosvd_function',
    'hosvd',
    'tprod',
    'ndim_unfold',
    'ndim_fold', 
    'svdtrunc',
    'cont_reconst',
    'dmd_reconst',
    'calculate_dmd_mode'
]