import numpy as np
import pennylane as qml
import sys,os
import torch
import torch.optim as optim
from torch.optim import lr_scheduler
from time import time
from pennylane_qchem.qchem import convert_observable 
from openfermion.utils import count_qubits
from transform import prepare_qubit_ham 
from clean_file import clean_all_except_py  

from pennylane.transforms.mitigate import fold_global, poly_extrapolate

import qiskit
import qiskit.providers.aer.noise as noise


# prepare the qubit hamilatonian
if (not os.path.exists('ham.npy') or not os.path.exists('fci.txt')):
#    clean_all_except_py()
    prepare_qubit_ham()
ham=np.load('ham.npy',allow_pickle=True).item()
H=convert_observable(ham)
n_qubits = count_qubits(ham)
fci=np.loadtxt('fci.txt')

# define the network structure
network=[n_qubits,1]

n_qubits_tot=sum(network)
ncycle=1

n_tot_params=6*network[0]
for icycle in range(ncycle):
    n_tot_params += 2*4*network[0]

print('The Architecture of network: ',network)
print('The cyle of each layer: ',ncycle)
print('The total parameters of network:',n_tot_params)

sys.stdout.flush()

#if os.path.exists('params.txt'):
#    theta=np.loadtxt('params.txt')
#else:
#    theta=np.random.uniform(high=2*np.pi,size=(n_tot_params))

def circuit(params):
   qml.BasisState(np.zeros(n_qubits_tot),wires=range(n_qubits_tot))

   nparams=0
   # for the input gates
   for i in range(network[0]):
      qml.RX(params[nparams+2*i], wires=i)
      qml.RY(params[nparams+2*i+1], wires=i)
      qml.RZ(params[nparams+2*i+2], wires=i)

   nparams=3*network[0]

   for icycle in range(ncycle):

      for i in range(network[0]):
         neuron_params = params[nparams:]
         # iterate over all input neurons and apply CAN gates
         for j in [ (i) % network[1] ]:
            tx1 = neuron_params[0]
            tx2 = neuron_params[1]
            ty1 = neuron_params[2]
            ty2 = neuron_params[3]
           
            qml.RX(tx1, wires=sum(network[:1])+j)
            qml.CNOT(wires=(sum(network[:1])+j,i))
            qml.RX(tx2, wires=i)

            qml.RZ(ty1, wires=sum(network[:1])+j)
            qml.CNOT(wires=(sum(network[:1])+j,i))
            qml.RZ(ty2, wires=i)

            qml.CNOT(wires=(sum(network[:1])+j,i))

         nparams += 4

      for i in range(network[0]):
         neuron_params = params[nparams:]
         # iterate over all input neurons and apply CAN gates
         for j in [ (i) % network[1] ]:
            tx1 = neuron_params[0]
            tx2 = neuron_params[1]
            ty1 = neuron_params[2]
            ty2 = neuron_params[3]
#            tz1 = neuron_params[4]
#            tz2 = neuron_params[5]

            qml.RX(tx1, wires=sum(network[:1])+j)
            qml.CNOT(wires=(sum(network[:1])+j,i))
            qml.RX(tx2, wires=i)

            qml.RZ(ty1, wires=sum(network[:1])+j)
            qml.CNOT(wires=(sum(network[:1])+j,i))
            qml.RZ(ty2, wires=i)

#            qml.RZ(2*tz1, wires=sum(network[:1])+j)
            qml.CNOT(wires=(sum(network[:1])+j,i))
#            qml.RZ(2*tz2, wires=i)
         nparams += 4

      for i in range(network[0]):
          qml.RX(params[nparams+2*i], wires=i)
          qml.RY(params[nparams+2*i+1], wires=i)
          qml.RZ(params[nparams+2*i+2], wires=i)
      nparams += 3*network[0]


dev=qml.device('default.qubit', wires=n_qubits_tot)
@qml.qnode(dev)
def cost_fn_no_noise(params):
    circuit(params)
    return qml.expval(H)


def generate_error_matrix(theta_,phi,delta_p,delta_s,delta_):
   error_mat=np.zeros((4,4),dtype=complex)
   error_mat[0][0]=1.0
   error_mat[1][1]=np.cos(theta_)*np.exp(complex(0+1j)*(delta_p+delta_s))
   error_mat[1][2]=complex(0-1j)*np.exp(complex(0+1j)*(delta_p-delta_))*np.sin(theta_)
   error_mat[2][1]=complex(0-1j)*np.exp(complex(0+1j)*(delta_p+delta_))*np.sin(theta_)
   error_mat[2][2]=np.cos(theta_)*np.exp(complex(0+1j)*(delta_p-delta_s))
   error_mat[3][3]=np.exp(complex(0+1j)*(2*delta_p-phi))

   return error_mat

nn=100

for i in range(nn):
   theta=np.random.uniform(high=2*np.pi,size=(n_tot_params))
   r_energy=cost_fn_no_noise(theta)

   theta_,phi,delta_p,delta_s,delta_=0.01*np.random.rand(5)
   error_matrix=generate_error_matrix(theta_,phi,delta_p,delta_s,delta_)
   
   error_1 = noise.mixed_unitary_error(noise_ops=[(error_matrix,1)])
   error_2 = noise.depolarizing_error(1-0.999, 1)
   error_3 = noise.depolarizing_error(1-0.994, 2)
   
   # Add errors to noise model
   noise_model = noise.NoiseModel()
   noise_model.add_all_qubit_quantum_error(error_1, ['cx'])
   noise_model.add_all_qubit_quantum_error(error_2, ['rx','rz'])
   noise_model.add_all_qubit_quantum_error(error_3, ['cx'])
   
   #print('noise_model:',noise_model)
   
   dev=qml.device('qiskit.aer', wires=n_qubits_tot, noise_model=noise_model,backend="aer_simulator_matrix_product_state")
   @qml.qnode(dev)
   def cost_fn(params):
      circuit(params)
      return qml.expval(H)
  
   r_energy_noise=cost_fn(theta) 

   with open('r_energy.txt','a') as f:
       np.savetxt(f,[[i,r_energy,r_energy_noise]])


#@qml.transforms.mitigate_with_zne([1., 2., 3.], fold_global, poly_extrapolate, extrapolate_kwargs={'order': 2})
##dev=qml.device('qiskit.aer', wires=n_qubits_tot, noise_model=noise_model,backend="aer_simulator_matrix_product_state")
#@qml.qnode(dev)
#def cost_with_zne(params):
#   circuit(params)
#   return qml.expval(H)
#
#print('final energy with noise(ZNE): ',cost_with_zne(theta))
#sys.stdout.flush()
