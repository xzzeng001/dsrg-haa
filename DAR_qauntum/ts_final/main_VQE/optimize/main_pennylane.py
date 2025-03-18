import numpy as np
import pennylane as qml
import sys,os
import torch
import torch.optim as optim
from torch.optim import lr_scheduler
from time import time
from pennylane_qchem.qchem import convert_observable 
from openfermion.utils import count_qubits

ham=np.load('ham.npy',allow_pickle=True).item()
H=convert_observable(ham)
n_qubits = count_qubits(ham)
#fci=np.loadtxt('fci.txt')

# define the network structure
network=[n_qubits,1]

n_qubits_tot=sum(network)
ncycle=2

n_tot_params=2*network[0]
for icycle in range(ncycle):
    n_tot_params += 4*network[0]

print('The Architecture of network: ',network)
print('The cyle of each layer: ',ncycle)
print('The total parameters of network:',n_tot_params)

sys.stdout.flush()

def circuit(params):
   qml.BasisState(np.zeros(n_qubits_tot),wires=range(n_qubits_tot))

   nparams=0
   # for the input gates
   for i in range(network[0]):
      qml.RX(params[nparams+2*i], wires=i)
#      qml.RY(params[nparams+2*i+1], wires=i)
      qml.RZ(params[nparams+2*i+2], wires=i)

   nparams=2*network[0]

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
            qml.RZ(ty1, wires=sum(network[:1])+j)

            qml.CNOT(wires=(sum(network[:1])+j,i))

            qml.RX(tx2, wires=i)
            qml.RZ(ty2, wires=i)

         nparams += 4


cpu_device = torch.device("cpu") 

dev = qml.device("default.qubit.torch", wires=n_qubits_tot)

@qml.qnode(dev,interface="torch", diff_method="backprop")
@qml.compile()
def cost_fn(params):
   circuit(params)
   return qml.expval(H)

#opt=qml.AdamOptimizer(stepsize=0.01)

if os.path.exists('params.txt'):
    theta=torch.tensor(np.loadtxt('params.txt'),requires_grad=True,dtype=torch.float64,device=cpu_device)
else:
    theta=torch.tensor(np.random.uniform(high=2*np.pi,size=(n_tot_params)),requires_grad=True,dtype=torch.float64,device=cpu_device)

max_iterations=10000
optimizer = torch.optim.Adam([theta],lr=0.01)
scheduler = lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.99)

prev_energy=cost_fn(theta).detach().numpy()
for n in range(max_iterations):
   r_energy=cost_fn(theta)
   energy=r_energy
   optimizer.zero_grad()

   energy.backward()

   optimizer.step()

   scheduler.step()

   r_energy=r_energy.detach().numpy()
#   conv=r_energy

   with open('energy.txt','a') as f:
       np.savetxt(f,[[n,r_energy]])

   if r_energy < prev_energy:
       np.savetxt("params.txt",[theta.detach().numpy()])
       prev_energy=r_energy

#if (not os.path.exists('1.png')):
#   qml.drawer.use_style('default')
#   fig, ax = qml.draw_mpl(cost_fn,decimals=4)(theta)
#   fig.savefig('1.png')

