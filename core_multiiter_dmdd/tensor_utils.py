import numpy as np
from typing import List, Tuple, Optional


def ndim_unfold(T: np.ndarray, dim: int) -> np.ndarray:
    """
    Layout an array in a given dimension to a matrix.
    
    Parameters:
    -----------
    T : np.ndarray
        Multidimensional array
    dim : int
        Dimension which becomes the column of the resulting matrix (0-indexed)
        
    Returns:
    --------
    np.ndarray
        Matrix with every vector of T along the dim dimension
    """
    ndim = T.ndim
    
    if ndim == 2:
        if dim == 1:
            return T.T
        else:
            return T
    else:
        # Move the specified dimension to the front
        axes = list(range(ndim))
        axes = [dim] + [i for i in axes if i != dim]
        T_permuted = np.transpose(T, axes)
        
        # Reshape to matrix
        siz = T_permuted.shape[0]
        M = T_permuted.reshape(siz, -1)
        return M


def ndim_fold(M: np.ndarray, dim: int, siz: List[int]) -> np.ndarray:
    """
    Restore a matrix into a multidimensional array with the given size.
    
    Parameters:
    -----------
    M : np.ndarray
        Matrix
    dim : int
        The columns of M will turn into this dimension of T (0-indexed)
    siz : List[int]
        Target tensor shape
        
    Returns:
    --------
    np.ndarray
        Multidimensional array containing the column vectors of M along dim dimension
    """
    ndim = len(siz)
    
    if ndim == 2:
        if dim == 1:
            return M.T
        else:
            return M
    else:
        # Create new size with dim moved to front
        new_size = [siz[dim]] + [siz[i] for i in range(ndim) if i != dim]
        new_size[0] = M.shape[0]
        
        # Reshape matrix to tensor
        T = M.reshape(new_size)
        
        # Permute back to original order
        axes = list(range(ndim))
        axes = [i+1 if i < dim else (0 if i == dim else i) for i in range(ndim)]
        T = np.transpose(T, axes)
        
        return T


def svdtrunc(A: np.ndarray, ns: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Truncated SVD decomposition.
    
    Parameters:
    -----------
    A : np.ndarray
        Matrix
    ns : int, optional
        Number of singular values to keep
        
    Returns:
    --------
    tuple
        U (left singular vectors), sv (singular values)
    """
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    sv = S
    
    if ns is not None and ns < len(sv):
        sv = sv[:ns]
        U = U[:, :ns]
    
    return U, sv


def tprod(S: np.ndarray, U: List[Optional[np.ndarray]]) -> np.ndarray:
    """
    Tensor product of an ndim-array and multiple matrices.
    
    Parameters:
    -----------
    S : np.ndarray
        Tensor (multidimensional array)
    U : List[Optional[np.ndarray]]
        List containing matrices for each dimension of S
        
    Returns:
    --------
    np.ndarray
        Result of the product
    """
    T = S.copy()
    siz = list(S.shape)
    
    for i, Ui in enumerate(U):
        if Ui is not None and i < len(siz):
            siz[i] = Ui.shape[0]
            H = ndim_unfold(T, i)
            T = ndim_fold(Ui @ H, i, siz)
    
    return T


def hosvd(T: np.ndarray, n: List[int]) -> Tuple[np.ndarray, np.ndarray, List[np.ndarray], List[np.ndarray], List[int]]:
    """
    High Order SVD of a multidimensional array.
    
    Parameters:
    -----------
    T : np.ndarray
        Multidimensional array
    n : List[int]
        Number of singular values to retain for each mode
        
    Returns:
    --------
    tuple
        TT (reconstructed tensor), S (core tensor), U (mode matrices), 
        sv (singular values), nn (adjusted dimensions)
    """
    M = list(T.shape)
    P = len(M)
    
    U = []
    sv = []
    nn = n.copy()
    
    # Adjust dimensions based on total size constraint
    produto = np.prod(nn)
    for i in range(P):
        nn[i] = min(nn[i], produto // nn[i] if nn[i] > 0 else M[i])
    
    for i in range(P):
        nn[i] = min(nn[i], produto // nn[i] if nn[i] > 0 else M[i])
        
        # Unfold tensor along dimension i
        A = ndim_unfold(T, i)
        
        # SVD-based reduction
        Ui, svi = svdtrunc(A, nn[i])
        
        if nn[i] < 2:
            U.append(np.column_stack([Ui[:, 0], np.zeros_like(Ui[:, 0])]))
        else:
            U.append(Ui[:, :nn[i]])
        
        sv.append(svi)
    
    # Compute core tensor
    UT = [Ui.T for Ui in U]
    S = tprod(T, UT)
    
    # Reconstruct tensor
    TT = tprod(S, U)
    
    return TT, S, U, sv, nn