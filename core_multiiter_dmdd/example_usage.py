"""
Example usage of the Core Multi-Iterative HODMD implementation.

This example demonstrates how to use the MultiIterativeHODMD class
for analyzing multi-dimensional tensor data.
"""

import numpy as np
from multi_iter_hodmd import MultiIterativeHODMD


def create_synthetic_data():
    """
    Create synthetic tensor data for testing.
    
    Returns:
    --------
    np.ndarray
        Synthetic 4D tensor (spatial_x, spatial_y, spatial_z, time)
    """
    # Create a synthetic 3D+time tensor
    nx, ny, nz, nt = 20, 15, 10, 100
    
    # Spatial grids
    x = np.linspace(0, 2*np.pi, nx)
    y = np.linspace(0, 2*np.pi, ny) 
    z = np.linspace(0, np.pi, nz)
    t = np.linspace(0, 10, nt)
    
    # Create meshgrids
    X, Y, Z, T = np.meshgrid(x, y, z, t, indexing='ij')
    
    # Create synthetic data with multiple modes
    # Mode 1: traveling wave
    mode1 = np.sin(X + 0.5*T) * np.cos(Y) * np.exp(-0.1*T)
    
    # Mode 2: oscillating pattern
    mode2 = 0.5 * np.cos(2*X) * np.sin(Y + T) * np.sin(Z) * np.exp(-0.05*T)
    
    # Mode 3: decaying mode
    mode3 = 0.3 * np.sin(X + Y + Z) * np.exp(-0.2*T) * np.cos(2*T)
    
    # Combine modes and add noise
    tensor_data = mode1 + mode2 + mode3 + 0.02 * np.random.randn(nx, ny, nz, nt)
    
    return tensor_data


def run_example():
    """
    Run the Multi-Iterative HODMD analysis example.
    """
    print("Creating synthetic data...")
    tensor_data = create_synthetic_data()
    print(f"Tensor shape: {tensor_data.shape}")
    
    # Initialize the HODMD solver
    print("\nInitializing Multi-Iterative HODMD solver...")
    hodmd = MultiIterativeHODMD(max_iterations=50)
    
    # Define parameters for analysis
    tolerances = [1e-3, 1e-2]  # Reduced set for faster testing
    d_values = [10, 20, 30]    # Reduced set for faster testing
    
    print("\nRunning Multi-Iterative HODMD analysis...")
    results = hodmd.fit(
        Tensor=tensor_data,
        time_pos=3,  # Last dimension is time (0-indexed)
        delta_t=0.1,
        tolerances=tolerances,
        d_values=d_values
    )
    
    # Display results summary
    print("\n" + "="*60)
    print("RESULTS SUMMARY")
    print("="*60)
    
    summary = hodmd.get_results_summary()
    for key, result in summary.items():
        print(f"\n{key}:")
        print(f"  Final iteration: {result['final_iteration']}")
        print(f"  RRMSE: {result['RRMSE']:.6e}")
        print(f"  Number of modes: {result['n_modes']}")
        print(f"  Parameters: d={result['parameters']['d']}, "
              f"ε₁={result['parameters']['varepsilon1']:.0e}, "
              f"ε₂={result['parameters']['varepsilon2']:.0e}")
    
    # Show detailed results for best case (lowest RRMSE)
    best_key = min(summary.keys(), key=lambda k: summary[k]['RRMSE'])
    best_result = results[best_key]
    
    print(f"\n" + "="*60)
    print(f"DETAILED RESULTS FOR BEST CASE: {best_key}")
    print("="*60)
    
    print(f"Growth rates: {best_result['GrowthRate']}")
    print(f"Frequencies: {best_result['Frequency']}")
    print(f"Amplitudes: {best_result['Amplitude']}")
    print(f"Final RRMSE: {best_result['RRMSE']:.6e}")
    print(f"DMD modes shape: {best_result['DMDmode'].shape}")
    print(f"Reconstructed tensor shape: {best_result['TensorReconst'].shape}")
    
    print(f"\nAnalysis completed successfully!")
    
    return results


if __name__ == "__main__":
    results = run_example()