import numpy as np
from typing import List
from tensor_utils import tprod


def cont_reconst(t: float, t0: float, u: np.ndarray, Time: np.ndarray, 
                 GrowthRate: np.ndarray, Frequency: np.ndarray) -> np.ndarray:
    """
    Continuous reconstruction using DMD expansion.
    
    Parameters:
    -----------
    t : float
        Current time
    t0 : float
        Initial time
    u : np.ndarray
        DMD modes
    Time : np.ndarray
        Time vector
    GrowthRate : np.ndarray
        Growth rates
    Frequency : np.ndarray
        Frequencies
        
    Returns:
    --------
    np.ndarray
        Reconstructed state at time t
    """
    N, M = u.shape
    vv = np.zeros(M, dtype=complex)
    
    for m in range(M):
        vv[m] = np.exp((GrowthRate[m] + 1j * Frequency[m]) * (t - t0))
    
    Reconst = u @ vv
    return Reconst


def dmd_reconst(GrowthRate: np.ndarray, Frequency: np.ndarray, hatMode: np.ndarray, 
                Time: np.ndarray, U: List[np.ndarray], S: np.ndarray, 
                sv: List[np.ndarray], nn: List[int], TimePos: int) -> np.ndarray:
    """
    Reconstruct the original tensor using the DMD expansion.
    
    Parameters:
    -----------
    GrowthRate : np.ndarray
        Growth rates of DMD modes
    Frequency : np.ndarray  
        Frequencies of DMD modes
    hatMode : np.ndarray
        Reduced DMD modes
    Time : np.ndarray
        Time vector
    U : List[np.ndarray]
        SVD modes from HOSVD
    S : np.ndarray
        Tensor core
    sv : List[np.ndarray]
        Singular values
    nn : List[int]
        Number of retained singular values
    TimePos : int
        Position of temporal dimension (0-indexed)
        
    Returns:
    --------
    np.ndarray
        Reconstructed tensor
    """
    N, _ = hatMode.shape
    K = len(Time)
    hatTReconst = np.zeros((N, K), dtype=complex)
    
    # Reconstruction using DMD expansion
    for k in range(K):
        hatTReconst[:, k] = cont_reconst(Time[k], Time[0], hatMode, Time, GrowthRate, Frequency)
    
    # Reconstruction of original tensor using reduced tensor and tensor core
    Unondim = U.copy()
    UTnondim_temp = np.zeros((nn[TimePos], hatTReconst.shape[1]), dtype=complex)
    
    for kk in range(nn[TimePos]):
        if sv[TimePos][kk] != 0:
            UTnondim_temp[kk, :] = hatTReconst[kk, :] / sv[TimePos][kk]
        else:
            UTnondim_temp[kk, :] = 0
    
    Unondim[TimePos] = UTnondim_temp.T
    
    TensorReconst = tprod(S, Unondim)
    
    return TensorReconst


def calculate_dmd_mode(N: int, hatMode: np.ndarray, Amplitude: np.ndarray, 
                       U: List[np.ndarray], S: np.ndarray, nn: List[int], 
                       TimePos: int) -> np.ndarray:
    """
    Calculate DMD modes in the original tensor space.
    
    Parameters:
    -----------
    N : int
        Temporal dimension size
    hatMode : np.ndarray
        Reduced DMD modes  
    Amplitude : np.ndarray
        Mode amplitudes
    U : List[np.ndarray]
        SVD modes from HOSVD
    S : np.ndarray
        Tensor core
    nn : List[int]
        Number of retained singular values
    TimePos : int
        Position of temporal dimension (0-indexed)
        
    Returns:
    --------
    np.ndarray
        DMD modes in original space
    """
    hatMode_m = np.zeros((N, len(Amplitude)))
    for ii in range(len(Amplitude)):
        if Amplitude[ii] != 0:
            hatMode_m[:, ii] = hatMode[:, ii] / Amplitude[ii]
        else:
            hatMode_m[:, ii] = 0
    
    ModesT = np.zeros((nn[TimePos], hatMode_m.shape[1]))
    for kk in range(nn[TimePos]):
        ModesT[kk, :] = hatMode_m[kk, :]
    
    # Temporal DMD modes in reduced dimension
    Modes = U.copy()
    Modes[TimePos] = ModesT.T
    
    # Reconstruction of temporal DMD modes
    DMDmode = tprod(S, Modes)
    
    return DMDmode