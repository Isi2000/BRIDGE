import numpy as np
from typing import List, Tuple
from tensor_utils import hosvd


def hosvd_function(Tensor: np.ndarray, varepsilon1: float, nn: List[int], 
                   n: List[int], TimePos: int) -> Tuple[np.ndarray, List[np.ndarray], np.ndarray, List[np.ndarray], List[int]]:
    """
    Perform HOSVD decomposition to calculate the reduced temporal matrix hatT.
    
    Parameters:
    -----------
    Tensor : np.ndarray
        Initial data tensor
    varepsilon1 : float
        Tolerance for SVD
    nn : List[int]
        Current number of singular values
    n : List[int]
        Initial tensor dimensions
    TimePos : int
        Position of temporal dimension (0-indexed)
        
    Returns:
    --------
    tuple
        hatT (reduced temporal matrix), U (temporal SVD modes), 
        S (tensor core), sv (singular values), nn1 (updated dimensions)
    """
    # Perform HOSVD and retain all singular values
    TT, S, U, sv, n_updated = hosvd(Tensor, n)
    
    # Set truncation of singular values using varepsilon1
    nn1 = nn.copy()
    
    for i in range(1, len(nn1)):
        count = 0
        if len(sv[i]) > 0:
            for j in range(len(sv[i])):
                if sv[i][j] / sv[i][0] >= varepsilon1:
                    count += 1
                else:
                    break
        nn1[i] = count
    
    # HOSVD retaining nn1 singular values
    TT, S, U, sv, nn1 = hosvd(Tensor, nn1)
    
    # Construct reduced temporal matrix
    UT = [Ui.T for Ui in U]
    
    hatT = np.zeros((nn1[TimePos], Tensor.shape[TimePos]))
    for kk in range(nn1[TimePos]):
        hatT[kk, :] = sv[TimePos][kk] * UT[TimePos][kk, :]
    
    return hatT, U, S, sv, nn1