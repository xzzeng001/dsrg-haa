import numpy as np
import scipy
import openfermion

ham=np.load('ham.npy',allow_pickle=True).item()
ham_matrix = openfermion.get_sparse_operator(ham)
e, eigvectors = scipy.sparse.linalg.eigsh(ham_matrix, k=1, which="SA")
print('fci: ',e)
