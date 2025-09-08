import numpy as np
from typing import Tuple


def dmdd_core(d: int, V: np.ndarray, time: np.ndarray, 
              varepsilon1: float, varepsilon: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    DMD with delay coordinates algorithm.
    
    Parameters:
    -----------
    d : int
        Parameter of DMD-d (higher order Koopman assumption)
    V : np.ndarray
        Snapshot matrix
    time : np.ndarray
        Time vector
    varepsilon1 : float
        First tolerance (SVD)
    varepsilon : float
        Second tolerance (DMD-d modes)
        
    Returns:
    --------
    tuple
        GrowthRate, Frequency, Amplitude, DMDmode
    """
    J, K = V.shape
    
    # SVD of the original data
    U, sigma, Vt = np.linalg.svd(V, full_matrices=False)
    sigmas = sigma
    n = len(sigmas)
    
    # Determine spatial complexity
    NormS = np.linalg.norm(sigmas)
    kk = 0
    for k in range(n):
        if np.linalg.norm(sigmas[k:]) / NormS > varepsilon1:
            kk += 1
    
    U = U[:, :kk]
    
    # Create reduced snapshots matrix
    hatT = np.diag(sigma[:kk]) @ Vt[:kk, :]
    N, _ = hatT.shape
    
    # Create the modified snapshot matrix
    tildeT = np.zeros((d * N, K - d + 1))
    for ppp in range(d):
        tildeT[ppp*N:(ppp+1)*N, :] = hatT[:, ppp:ppp+K-d+1]
    
    # Dimension reduction
    U1, sigma1, V1t = np.linalg.svd(tildeT, full_matrices=False)
    sigmas1 = sigma1
    
    Deltat = time[1] - time[0]
    n = len(sigmas1)
    
    NormS = np.linalg.norm(sigmas1)
    kk1 = 0
    for k in range(n):
        if np.linalg.norm(sigmas1[k:]) / NormS > varepsilon1:
            kk1 += 1
    
    U1 = U1[:, :kk1]
    hatT1 = np.diag(sigma1[:kk1]) @ V1t[:kk1, :]
    
    # Reduced modified snapshot matrix
    _, K1 = hatT1.shape
    tildeU1, tildeSigma, tildeU2t = np.linalg.svd(hatT1[:, :K1-1], full_matrices=False)
    
    # Reduced modified Koopman matrix
    tildeR = hatT1[:, 1:K1] @ tildeU2t.T @ np.linalg.inv(np.diag(tildeSigma)) @ tildeU1.T
    eigenvalues, tildeQ = np.linalg.eig(tildeR)
    
    M = len(eigenvalues)
    qq = np.log(eigenvalues)
    GrowthRate = np.real(qq) / Deltat
    Frequency = np.imag(qq) / Deltat
    
    Q = U1 @ tildeQ
    Q = Q[(d-1)*N:d*N, :]
    NN, MMM = Q.shape
    
    # Normalize Q
    for m in range(MMM):
        Q[:, m] = Q[:, m] / np.linalg.norm(Q[:, m])
    
    # Calculate amplitudes
    Mm = np.zeros((NN * K, M))
    Bb = np.zeros((NN * K, 1))
    aa = np.eye(MMM)
    tildeMM = np.diag(eigenvalues)
    
    for k in range(K):
        Mm[k*NN:(k+1)*NN, :] = Q @ aa
        aa = aa @ tildeMM
        Bb[k*NN:(k+1)*NN, 0] = hatT[:, k]
    
    Ur, Sigmar, Vr = np.linalg.svd(Mm, full_matrices=False)
    a = Vr.T @ np.linalg.inv(np.diag(Sigmar)) @ Ur.T @ Bb
    a = a.flatten()
    
    u = np.zeros((NN, M))
    for m in range(M):
        u[:, m] = a[m] * Q[:, m]
    
    Amplitude = np.zeros(M)
    for m in range(M):
        aca = U @ u[:, m]
        Amplitude[m] = np.linalg.norm(aca) / np.sqrt(J)
    
    # Sort by amplitude
    UU = np.vstack([u, GrowthRate.reshape(1, -1), Frequency.reshape(1, -1), Amplitude.reshape(1, -1)])
    sort_idx = np.argsort(-Amplitude)
    UU = UU[:, sort_idx]
    
    u = UU[:NN, :]
    GrowthRate = UU[NN, :]
    Frequency = UU[NN+1, :]
    Amplitude = UU[NN+2, :]
    
    # Determine spectral complexity
    kk3 = 0
    for m in range(M):
        if Amplitude[m] / Amplitude[0] > varepsilon:
            kk3 += 1
    
    u = u[:, :kk3]
    GrowthRate = GrowthRate[:kk3]
    Frequency = Frequency[:kk3]
    Amplitude = Amplitude[:kk3]
    
    # Calculate DMD modes
    DMDmode = np.zeros((J, kk3))
    for m in range(kk3):
        NormMode = np.linalg.norm(U @ u[:, m]) / np.sqrt(J)
        DMDmode[:, m] = U @ u[:, m] / NormMode
    
    return GrowthRate, Frequency, Amplitude, DMDmode