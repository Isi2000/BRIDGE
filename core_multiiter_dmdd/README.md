# Core Multi-Iterative Dynamic Mode Decomposition with Delay Coordinates

This Python package implements the Multi-dimensional Higher-Order Dynamic Mode Decomposition (HODMD) with iterative refinement, specifically designed for complex signal analysis without plotting functionality.

## Features

- **Multi-Iterative HODMD**: Iterative refinement algorithm for improved accuracy
- **DMD with Delay Coordinates**: Higher-order Koopman operator approximation
- **Tensor Operations**: Complete set of tensor manipulation utilities
- **HOSVD Decomposition**: Higher-Order Singular Value Decomposition
- **Reconstruction**: Full tensor reconstruction capabilities
- **No Plotting Dependencies**: Pure computational implementation

## Files Structure

```
core_multiiter_dmdd/
├── __init__.py                 # Package initialization
├── multi_iter_hodmd.py         # Main iterative algorithm class
├── dmdd_core.py               # Core DMD with delay coordinates
├── hosvd_function.py          # HOSVD wrapper function
├── tensor_utils.py            # Tensor manipulation utilities
├── reconstruction.py          # Reconstruction functions
├── example_usage.py           # Usage example
└── README.md                  # This file
```

## Core Components

### 1. MultiIterativeHODMD Class
Main class for performing multi-iterative HODMD analysis:
```python
from multi_iter_hodmd import MultiIterativeHODMD

hodmd = MultiIterativeHODMD(max_iterations=1000)
results = hodmd.fit(tensor_data, time_pos=3, delta_t=0.1)
```

### 2. DMD with Delay Coordinates
Core algorithm implementing DMD-d:
```python
from dmdd_core import dmdd_core

growth_rate, frequency, amplitude, modes = dmdd_core(d, V, time, eps1, eps2)
```

### 3. Tensor Utilities
Essential tensor operations:
```python
from tensor_utils import hosvd, tprod, ndim_unfold, ndim_fold

# Higher-order SVD
TT, S, U, sv, nn = hosvd(tensor, dimensions)

# Tensor product
result = tprod(core_tensor, mode_matrices)
```

### 4. Reconstruction Functions
For rebuilding tensors from DMD expansion:
```python
from reconstruction import dmd_reconst, calculate_dmd_mode

tensor_reconstructed = dmd_reconst(growth_rate, frequency, modes, time, U, S, sv, nn, time_pos)
```

## Algorithm Overview

The Multi-Iterative HODMD algorithm performs the following steps:

1. **HOSVD Decomposition**: Decomposes the input tensor using Higher-Order SVD
2. **DMD-d Analysis**: Applies DMD with delay coordinates to the reduced temporal matrix
3. **Mode Selection**: Filters modes based on amplitude thresholds
4. **Reconstruction**: Rebuilds the tensor using the DMD expansion
5. **Iteration**: Repeats the process using the reconstructed tensor until convergence

## Parameters

### Key Parameters:
- `d`: DMD delay parameter (higher-order Koopman assumption)
- `varepsilon1`: SVD tolerance for dimension reduction
- `varepsilon2`: Mode selection tolerance (amplitude threshold)
- `time_pos`: Position of temporal dimension in tensor (0-indexed)
- `delta_t`: Time step between snapshots

### Typical Values:
- `d`: 10-70 (depending on system complexity)
- `varepsilon1`: 1e-4 to 1e-2
- `varepsilon2`: 1e-4 to 1e-2
- `max_iterations`: 50-1000

## Usage Example

```python
import numpy as np
from multi_iter_hodmd import MultiIterativeHODMD

# Create or load your tensor data
# tensor_data.shape = (spatial_dim1, spatial_dim2, spatial_dim3, time)

# Initialize solver
hodmd = MultiIterativeHODMD(max_iterations=100)

# Run analysis
results = hodmd.fit(
    Tensor=tensor_data,
    time_pos=3,  # Time is last dimension
    delta_t=0.1,
    tolerances=[1e-4, 1e-3, 1e-2],
    d_values=[10, 20, 30, 40, 50]
)

# Get summary
summary = hodmd.get_results_summary()
print(summary)

# Access specific results
best_key = min(summary.keys(), key=lambda k: summary[k]['RRMSE'])
best_result = results[best_key]

print(f"Growth rates: {best_result['GrowthRate']}")
print(f"Frequencies: {best_result['Frequency']}")
print(f"Amplitudes: {best_result['Amplitude']}")
```

## Dependencies

- `numpy`: Numerical computations
- `scipy`: Scientific computing (SVD, eigendecomposition)
- `typing`: Type hints

## Output

The algorithm returns a comprehensive results dictionary containing:

- `GrowthRate`: Growth rates of identified modes
- `Frequency`: Frequencies of identified modes  
- `Amplitude`: Amplitudes of identified modes
- `DMDmode`: DMD modes in original tensor space
- `TensorReconst`: Reconstructed tensor
- `RRMSE`: Relative root mean square error
- `iteration_history`: Convergence history
- `parameters`: Analysis parameters used

## Mathematical Background

This implementation is based on the Higher-Order Dynamic Mode Decomposition algorithm presented in:

*Le Clainche, S., Vega, J.M. & Soria, J., "Higher order dynamic mode decomposition of noisy experimental data: The flow structure of a zero-net-mass-flux jet", Experimental Thermal and Fluid Science, 2018, 88, pp. 336-353*

And the foundational HODMD work:

*Le Clainche, S. & Vega, J.M., "Higher order dynamic mode decomposition", SIAM J. Appl. Dyn. Sys., 16(2), pp. 882-925, 2017*

## Performance Notes

- Memory usage scales with tensor size and number of modes
- Computational complexity depends on tensor dimensions and iteration count
- For large tensors, consider reducing tolerance values or limiting d_values range
- The algorithm automatically handles convergence detection

## Differences from MATLAB Version

This Python implementation:
- Removes all plotting functionality
- Uses 0-based indexing (vs 1-based in MATLAB)
- Provides object-oriented interface
- Includes comprehensive type hints
- Returns structured results dictionary
- Handles errors gracefully with NumPy/SciPy