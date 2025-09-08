import numpy as np
from typing import List, Tuple, Optional
from dmdd_core import dmdd_core
from hosvd_function import hosvd_function
from reconstruction import dmd_reconst, calculate_dmd_mode


class MultiIterativeHODMD:
    """
    Multi-dimensional Higher-Order Dynamic Mode Decomposition with iterative refinement.
    """
    
    def __init__(self, max_iterations: int = 1000):
        """
        Initialize the Multi-Iterative HODMD solver.
        
        Parameters:
        -----------
        max_iterations : int
            Maximum number of iterations
        """
        self.max_iterations = max_iterations
        self.results = {}
    
    def fit(self, Tensor: np.ndarray, time_pos: int = -1, delta_t: float = 1.0, 
            tolerances: List[float] = None, d_values: List[int] = None) -> dict:
        """
        Fit the Multi-Iterative HODMD model to the tensor data.
        
        Parameters:
        -----------
        Tensor : np.ndarray
            Input tensor data
        time_pos : int
            Position of time dimension (default: last dimension)
        delta_t : float
            Time step between snapshots
        tolerances : List[float]
            Tolerance values to test (default: [1e-4, 1e-3, 1e-2])
        d_values : List[int]
            DMD delay values to test (default: [10, 20, 30, 40, 50, 60, 70])
            
        Returns:
        --------
        dict
            Results dictionary containing all computed solutions
        """
        if tolerances is None:
            tolerances = [1e-4, 1e-3, 1e-2]
        if d_values is None:
            d_values = [10, 20, 30, 40, 50, 60, 70]
        
        # Handle time position
        if time_pos == -1:
            time_pos = Tensor.ndim - 1
        
        # Create time vector
        snap_count = Tensor.shape[time_pos]
        time = np.arange(1, snap_count + 1) * delta_t
        
        results = {}
        
        for tol_idx, tolerance in enumerate(tolerances):
            for d_idx, d in enumerate(d_values):
                print(f"Processing tolerance: {tolerance}, d: {d}")
                
                result = self._run_iterative_algorithm(
                    Tensor, time, time_pos, d, tolerance, tolerance
                )
                
                key = f"tol_{tolerance:.0e}_d_{d}"
                results[key] = result
        
        self.results = results
        return results
    
    def _run_iterative_algorithm(self, Tensor: np.ndarray, time: np.ndarray, 
                                 time_pos: int, d: int, varepsilon1: float, 
                                 varepsilon2: float) -> dict:
        """
        Run the iterative HODMD algorithm.
        
        Parameters:
        -----------
        Tensor : np.ndarray
            Input tensor
        time : np.ndarray
            Time vector
        time_pos : int
            Position of time dimension
        d : int
            DMD delay parameter
        varepsilon1 : float
            First tolerance (SVD)
        varepsilon2 : float
            Second tolerance (DMD modes)
            
        Returns:
        --------
        dict
            Results for this parameter combination
        """
        Tensor0 = Tensor.copy()
        
        # Initialize dimensions
        nn0 = list(Tensor.shape)
        nn = [nn0[0]] + [0] * (len(nn0) - 1)
        
        # Store results
        iteration_results = []
        
        for iteration in range(self.max_iterations):
            print(f"Iteration: {iteration + 1}")
            
            if iteration > 0:
                Tensor = TensorReconst
            
            # Perform HOSVD decomposition
            hatT, U, S, sv, nn1 = hosvd_function(Tensor, varepsilon1, nn, nn0, time_pos)
            
            # Perform DMD-d on reduced temporal matrix
            if d > 1:
                GrowthRate, Frequency, Amplitude, hatMode = dmdd_core(
                    d, hatT, time, varepsilon1, varepsilon2
                )
            else:
                # For d=1, use standard DMD (simplified version)
                GrowthRate, Frequency, Amplitude, hatMode = self._standard_dmd(
                    hatT, time, varepsilon1, varepsilon2
                )
            
            # Reconstruct original tensor using DMD expansion
            TensorReconst = dmd_reconst(GrowthRate, Frequency, hatMode, time, U, S, sv, nn1, time_pos)
            TensorReconst = np.real(TensorReconst)
            
            # Calculate relative mean square error
            RRMSE = np.linalg.norm(Tensor.flatten() - TensorReconst.flatten()) / np.linalg.norm(Tensor.flatten())
            
            # Store iteration results
            iter_result = {
                'iteration': iteration + 1,
                'GrowthRate': GrowthRate.copy(),
                'Frequency': Frequency.copy(),
                'Amplitude': Amplitude.copy(),
                'RRMSE': RRMSE,
                'nn': nn1.copy()
            }
            iteration_results.append(iter_result)
            
            print(f"RRMSE: {RRMSE:.6e}")
            print(f"Modes found: {len(GrowthRate)}")
            
            # Check convergence (number of singular values unchanged)
            if len(nn) == len(nn1):
                converged = all(nn[i] == nn1[i] for i in range(1, len(nn1)))
                if converged:
                    print(f"Converged after {iteration + 1} iterations")
                    break
            
            nn = nn1.copy()
        
        # Calculate DMD modes
        N, _ = hatT.shape
        DMDmode = calculate_dmd_mode(N, hatMode, Amplitude, U, S, nn1, time_pos)
        
        final_result = {
            'final_iteration': iteration + 1,
            'GrowthRate': GrowthRate,
            'Frequency': Frequency, 
            'Amplitude': Amplitude,
            'DMDmode': DMDmode,
            'TensorReconst': TensorReconst,
            'RRMSE': RRMSE,
            'iteration_history': iteration_results,
            'parameters': {
                'd': d,
                'varepsilon1': varepsilon1,
                'varepsilon2': varepsilon2,
                'time_pos': time_pos
            }
        }
        
        return final_result
    
    def _standard_dmd(self, V: np.ndarray, time: np.ndarray, 
                      varepsilon1: float, varepsilon2: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Standard DMD algorithm for d=1 case.
        
        Parameters:
        -----------
        V : np.ndarray
            Snapshot matrix
        time : np.ndarray
            Time vector
        varepsilon1 : float
            SVD tolerance
        varepsilon2 : float
            Mode selection tolerance
            
        Returns:
        --------
        tuple
            GrowthRate, Frequency, Amplitude, DMDmode
        """
        # SVD of snapshot matrix
        U, sigma, Vt = np.linalg.svd(V, full_matrices=False)
        
        # Truncate based on tolerance
        sigmas = sigma
        NormS = np.linalg.norm(sigmas)
        kk = 0
        for k in range(len(sigmas)):
            if np.linalg.norm(sigmas[k:]) / NormS > varepsilon1:
                kk += 1
        
        U = U[:, :kk]
        sigma = sigma[:kk]
        Vt = Vt[:kk, :]
        
        # Build Koopman matrix
        X1 = U.T @ V[:, :-1]
        X2 = U.T @ V[:, 1:]
        A_tilde = X2 @ np.linalg.pinv(X1)
        
        # Eigendecomposition
        eigenvalues, W = np.linalg.eig(A_tilde)
        
        # Calculate growth rates and frequencies
        dt = time[1] - time[0]
        GrowthRate = np.real(np.log(eigenvalues)) / dt
        Frequency = np.imag(np.log(eigenvalues)) / dt
        
        # Calculate modes and amplitudes
        modes = U @ W
        
        # Calculate amplitudes using least squares
        b = V[:, 0]
        amplitudes = np.linalg.lstsq(modes, b, rcond=None)[0]
        Amplitude = np.abs(amplitudes)
        
        # Sort by amplitude
        sort_idx = np.argsort(-Amplitude)
        GrowthRate = GrowthRate[sort_idx]
        Frequency = Frequency[sort_idx]
        Amplitude = Amplitude[sort_idx]
        modes = modes[:, sort_idx]
        
        # Select modes based on amplitude threshold
        n_modes = np.sum(Amplitude / Amplitude[0] > varepsilon2)
        
        return (GrowthRate[:n_modes], Frequency[:n_modes], 
                Amplitude[:n_modes], modes[:, :n_modes])
    
    def get_results_summary(self) -> dict:
        """
        Get a summary of all computed results.
        
        Returns:
        --------
        dict
            Summary of results across all parameter combinations
        """
        if not self.results:
            return {}
        
        summary = {}
        for key, result in self.results.items():
            summary[key] = {
                'final_iteration': result['final_iteration'],
                'RRMSE': result['RRMSE'],
                'n_modes': len(result['GrowthRate']),
                'parameters': result['parameters']
            }
        
        return summary