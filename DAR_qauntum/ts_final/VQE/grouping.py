import numpy as np
#import openfermion
#import scipy
#from openfermion.utils import count_qubits
from pennylane_qchem.qchem import convert_observable
import pennylane as qml
import os

from typing import Union, Optional, List, Tuple, Any, Dict

def ham_grouping(ham: Any):
    '''
    Grouping the hamiltonian based on the anticommutation relation
    mainly based on the penneylane
    '''

    os.mkdir('group')
    file_ham='group/ham.txt'
    file_ham_groups='group/ham_grouping.txt'

#    n_qubits = count_qubits(ham)
    nterms = len(ham.terms)
    H=convert_observable(ham)

    coeffs=H.coeffs
    obs=H.ops
    with open(file_ham,'ab') as f:
        np.savetxt(f,coeffs)
        np.savetxt(f,obs,fmt='%s')

    obs_groupings, coeffs_groupings=qml.grouping.group_observables(obs, coeffs, 'anticommuting', 'lf')
    with open(file_ham_groups,'ab') as f:
        np.savetxt(f,coeffs_groupings,fmt='%s')
        np.savetxt(f,obs_groupings,fmt='%s')

#    print('obs_groupings: ',coeffs_groupings)
#    print('coeffs_groupings.shape: ',len(coeffs_groupings))
#    ham_matrix = openfermion.get_sparse_operator(ham)
#    e, eigvectors = scipy.sparse.linalg.eigsh(ham_matrix, k=1, which="SA")
#
#    with open(outfile,'a') as f:
#        np.savetxt(f,[[ii,n_qubits,nterms,len(coeffs_groupings),e]])

#    print('ii, n_qubits, nterms, ngroup, e:',ii,n_qubits,nterms,len(coeffs_groupings),e)
