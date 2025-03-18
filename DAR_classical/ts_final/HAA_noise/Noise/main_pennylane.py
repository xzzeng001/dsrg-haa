#from pennylane import qchem
import numpy as np
import pennylane as qml
import sys,os
#import torch
#import torch.optim as optim
#from torch.optim import lr_scheduler
import pyscf
import pyscf.cc
import pyscf.fci
from time import time
from pyscf_func import obtain_Hamiltonian
from pennylane_qchem.qchem import convert_observable 
from openfermion.utils import count_qubits

from functools import partial
from pennylane.transforms import fold_global, poly_extrapolate

import qiskit
import qiskit_aer.noise as noise


if os.path.exists('ham.npy'):
    ham=np.load('ham.npy',allow_pickle=True).item()
    H=convert_observable(ham)
    n_qubits = count_qubits(ham)
else:
    H=convert_observable(ham)
    np.save('ham.npy',ham)

network=[n_qubits,2]

n_qubits_tot=sum(network)
ncycle=3

n_tot_params=3*network[0]
for icycle in range(ncycle):
    n_tot_params += 3*network[0]*network[1]

print('The Architecture of network: ',network)
print('The cyle of each layer: ',ncycle)
print('The total parameters of network:',n_tot_params)

sys.stdout.flush()

if os.path.exists('params.txt'):
    theta=np.loadtxt('params.txt')
else:
    theta=np.random.uniform(high=2*np.pi,size=(n_tot_params))

def circuit(params):
   qml.BasisState(np.zeros(n_qubits_tot),wires=range(n_qubits_tot))

   nparams=0
   # for the input gates
   for i in range(network[0]):
      qml.U3(params[nparams+i*3], params[nparams+i*3+1], params[nparams+i*3+2], wires=i)

   nparams=3*network[0]

   for icycle in range(ncycle):

      # for the intermediate layer
      for i in range(network[0]):
         # parameters of the respective "neuron gates"
         # (can be larer than needed, overflow will be ignored)
         neuron_params = params[nparams:]
         # iterate over all input neurons and apply CAN gates
         for j in range(network[1]):
            tx, ty, tz = neuron_params[j*3:(j+1)*3]
            qml.IsingXX(2*tx, wires=(sum(network[:1])+j,i))
            qml.IsingYY(2*ty, wires=(sum(network[:1])+j,i))
            qml.IsingZZ(2*tz, wires=(sum(network[:1])+j,i))
#            print('two target:',sum(network[:n_layer-ilayer-2])+j,sum(network[:n_layer-ilayer-1])+i)
         nparams += 3*network[1]


dev=qml.device('default.qubit', wires=n_qubits_tot)
@qml.qnode(dev)
def cost_fn_no_noise(params):
    circuit(params)
    return qml.expval(H)

e_no_noise=cost_fn_no_noise(theta)
print('final energy(without noise): ',e_no_noise)
sys.stdout.flush()

with open('xx_energy.txt','a') as f:
    np.savetxt(f,[[0,e_no_noise,e_no_noise]])


def generate_error_matrix(theta_,phi,delta_p,delta_s,delta_):
   error_mat=np.zeros((4,4),dtype=complex)
   error_mat[0][0]=1.0
   error_mat[1][1]=np.cos(theta_)*np.exp(complex(0+1j)*(delta_p+delta_s))
   error_mat[1][2]=complex(0-1j)*np.exp(complex(0+1j)*(delta_p-delta_))*np.sin(theta_)
   error_mat[2][1]=complex(0-1j)*np.exp(complex(0+1j)*(delta_p+delta_))*np.sin(theta_)
   error_mat[2][2]=np.cos(theta_)*np.exp(complex(0+1j)*(delta_p-delta_s)) 
   error_mat[3][3]=np.exp(complex(0+1j)*(2*delta_p-phi))

   return error_mat


theta_,phi,delta_p,delta_s,delta_=0.01*np.random.rand(5)
error_matrix=generate_error_matrix(theta_,phi,delta_p,delta_s,delta_)

error_1 = noise.mixed_unitary_error(noise_ops=[(error_matrix,1)])
error_2 = noise.depolarizing_error(1-0.999, 1)

for xx in np.arange(0.9,0.9999,0.002):
    error_3 = noise.depolarizing_error(1-xx, 2)

    # Add errors to noise model
    noise_model = noise.NoiseModel()
#    noise_model.add_all_qubit_quantum_error(error_1, ['cx'])
    noise_model.add_all_qubit_quantum_error(error_2, ['u3'])
    noise_model.add_all_qubit_quantum_error(error_3, ['rxx','ryy','rzz'])
    
    dev1=qml.device('qiskit.aer',wires=n_qubits_tot,noise_model=noise_model)
    @qml.qnode(dev1)
    def cost_fn(params):
       circuit(params)
       return qml.expval(H)
    
    e_noise=cost_fn(theta)
    print('final energy with noise: ',e_noise)
    sys.stdout.flush()
    
    
    #dev=qml.transforms.insert()(dev1)

#    @partial(qml.transforms.mitigate_with_zne,[1., 2., 3.], fold_global, poly_extrapolate, extrapolate_kwargs={'order': 1})
#    #dev=qml.device('qiskit.aer', wires=n_qubits_tot, noise_model=noise_model)
#    @qml.qnode(dev)
#    def cost_with_zne(params):
#       circuit(params)
#       return qml.expval(H)
#    
#    e_noise_zne=cost_with_zne(theta)
#    print('final energy with noise(ZNE): ',e_noise_zne)
#    sys.stdout.flush()

    with open('xx_energy.txt','a') as f:
        np.savetxt(f,[[xx,e_noise]])#,e_noise_zne]])

