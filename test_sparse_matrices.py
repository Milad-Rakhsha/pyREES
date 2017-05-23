import numpy as np
import scipy.sparse as sparse
import REESMath.matrix3 as M

if __name__ == '__main__':

    # rows = np.array([0, 1, 2, 3], dtype=np.int32)
    # cols = np.array([0, 1, 2, 3], dtype=np.int32)
    # vals = np.array([0, 1, 2, 3], dtype=np.float64)
    #
    # A = sparse.coo_matrix((vals, (rows, cols)), shape=(4, 4))
    # B = A.tocsr()
    #
    # print(A, type(A))
    # print(B, type(B))
    #
    # D1 = M.diag(1,1,1)
    # D2 = M.diag(2,2,2)
    # D3 = M.diag(3,3,3)
    #
    # diagonal_blocks = [D1, D2, D3]
    #
    # D = sparse.block_diag(diagonal_blocks, format='csr')
    # print(D, type(D))
    #
    # K = D.tobsr(blocksize=(3,3))
    #
    # print(K, type(K))
