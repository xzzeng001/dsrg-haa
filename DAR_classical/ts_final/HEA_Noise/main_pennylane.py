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

# prepare the qubit hamilatonian
#if (not os.path.exists('ham.npy') or not os.path.exists('fci.txt')):
#    clean_all_except_py()
#    prepare_qubit_ham()

ham=np.load('ham.npy',allow_pickle=True).item()
H=convert_observable(ham)
n_qubits = count_qubits(ham)
fci=np.loadtxt('fci.txt')

# define the network structure
network=[n_qubits,0]

n_qubits_tot=sum(network)
ncycle=4

n_tot_params=3*network[0]+3*ncycle*(network[0]-1)
#for icycle in range(ncycle):
#    n_tot_params += 3*network[0]*network[1]

print('The Architecture of network: ',network)
print('The cyle of each layer: ',ncycle)
print('The total parameters of network:',n_tot_params)

sys.stdout.flush()

def circuit(params):
   qml.BasisState(np.zeros(n_qubits_tot),wires=range(n_qubits_tot))

   nparams=0
   # for the input gates
   for i in range(network[0]):
      qml.U3(params[nparams+i*3], params[nparams+i*3+1], params[nparams+i*3+2], wires=i)

   nparams=3*network[0]

   # for the intermediate layer
   for icycle in range(ncycle):
      for i in range(network[0]-1):
         neuron_params = params[nparams:]
         tx, ty, tz = neuron_params[0:3]
         qml.IsingXX(2*tx, wires=(i,i+1))
         qml.IsingYY(2*ty, wires=(i,i+1))
         qml.IsingZZ(2*tz, wires=(i,i+1))
         nparams += 3

#assert torch.cuda.is_available()
cpu_device = torch.device("cpu") 

dev = qml.device("default.qubit.torch", wires=n_qubits_tot)

@qml.qnode(dev,interface="torch", diff_method="backprop")
def cost_fn(params):
   circuit(params)
   return qml.expval(H)

if os.path.exists('params.txt'):
    theta=torch.tensor(np.loadtxt('params.txt'),requires_grad=True,dtype=torch.float64,device=cpu_device)
else:
    theta=torch.tensor(np.random.uniform(high=2*np.pi,size=(n_tot_params)),requires_grad=True,dtype=torch.float64,device=cpu_device)

max_iterations=10000
optimizer = torch.optim.Adam([theta],lr=0.01)
scheduler = lr_scheduler.StepLR(optimizer, step_size=100, gamma=0.99)

prev_energy=cost_fn(theta).detach().numpy()
for n in range(max_iterations):
#   theta,prev_energy=opt.step_and_cost(cost_fn,theta)
   r_energy=cost_fn(theta)
   energy=r_energy
   optimizer.zero_grad()
   energy.backward()
   optimizer.step()

   scheduler.step()

   r_energy=r_energy.detach().numpy()
   conv=r_energy

   with open('energy.txt','a') as f:
       np.savetxt(f,[[n,r_energy,r_energy-fci]])

   if r_energy < prev_energy:
       np.savetxt("params.txt",[theta.detach().numpy()])
       prev_energy=r_energy

#if (not os.path.exists('1.png')):
#   qml.drawer.use_style('default')
#   fig, ax = qml.draw_mpl(cost_fn,decimals=4)(theta)
#   fig.savefig('1.png')

