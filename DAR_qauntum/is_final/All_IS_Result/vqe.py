import numpy as np
import pennylane as qml
from pennylane_qchem.qchem import convert_observable
from openfermion.utils import count_qubits
import os

def main_vqe(nn,params):

   # do the main vqe to calculate the energy

   # read the ham
   if os.path.exists('ham.npy'):
       ham=np.load('ham.npy',allow_pickle=True).item()
       H=convert_observable(ham)
       qubits = count_qubits(ham)

#   print('H: ',H)
   n_qubits_tot=qubits+nn

   # define the circuit:
   def circuit(params):
      qml.BasisState(np.zeros(n_qubits_tot),wires=range(n_qubits_tot))

      qml.Hadamard(wires=0)
      qml.RZ(params[0],wires=0)
      qml.Hadamard(wires=0)
      qml.RZ(params[1],wires=0)
      qml.Hadamard(wires=0)

      qml.Hadamard(wires=1)
      qml.RZ(params[2],wires=1)
      qml.Hadamard(1)
      qml.RZ(params[3],wires=1)
      qml.Hadamard(wires=1) 

      qml.Hadamard(wires=2)
      qml.RZ(params[4],wires=2)
      qml.Hadamard(wires=2)
      qml.RZ(params[5],wires=2)

      qml.CZ(wires=(0,2))

      qml.RZ(params[6],wires=0)
      qml.Hadamard(wires=0)
      qml.RZ(params[7],wires=0) 
      qml.Hadamard(wires=0)

      qml.Hadamard(wires=2)
      qml.RZ(params[8],wires=2)
      qml.Hadamard(wires=2)
      qml.RZ(params[9],wires=2)

      qml.CZ(wires=(1,2))

      qml.RZ(params[10],wires=1)
      qml.Hadamard(wires=1)
      qml.RZ(params[11],wires=1)
      qml.Hadamard(wires=1)

      qml.Hadamard(wires=2)
      qml.RZ(params[12],wires=1)
      qml.Hadamard(wires=2)
      qml.RZ(params[13],wires=1)

      qml.CZ(wires=(0,2))

      qml.RZ(params[14],wires=0)
      qml.Hadamard(wires=0)

      qml.Hadamard(wires=2)
      qml.RZ(params[15],wires=2)
      qml.Hadamard(wires=2)
      qml.RZ(params[16],wires=2)

      qml.CZ(wires=(1,2))

      qml.RZ(params[17],wires=1)
      qml.Hadamard(wires=1)

   
   
   dev = qml.device("default.qubit", wires=n_qubits_tot)
   
   @qml.qnode(dev)
   def cost_fn(params):
      circuit(params)
      return qml.expval(H)
  
   return cost_fn(params) 
