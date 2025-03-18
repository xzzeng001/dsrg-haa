from qiskit import QuantumCircuit, QuantumRegister, transpile, Aer, execute
from typing import Union, Optional, List, Tuple, Any, Dict
import numpy as np

# add one qubit gate
def add_one_qubit_gates(circ: QuantumCircuit,
                        q_reg: QuantumRegister,
                        network: Any,
                        params: Any,
                        iparams: Any) -> Tuple[QuantumCircuit,int]:

   for i in range(network[0]):
      circ.rx(params[iparams], q_reg[i])
      circ.ry(params[iparams+1], q_reg[i])
      circ.rx(params[iparams+2], q_reg[i])
      iparams += 3
      
   return circ,iparams

# add two-qubit gate
def add_two_qubit_gates(circ: QuantumCircuit,
                        q_reg: QuantumRegister,
                        network: Any,
                        params: Any,
                        iparams: Any) -> Tuple[QuantumCircuit,int]:

   for i in range(network[0]):
      for j in [(i) % network[1]]:
         circ.rx(params[iparams], q_reg[sum(network[:1])+j])
         circ.cnot(sum(network[:1])+j,i)
         circ.rx(params[iparams+1], q_reg[i]) 

         circ.rz(params[iparams+2], q_reg[sum(network[:1])+j])
         circ.cnot(sum(network[:1])+j,i)
         circ.rz(params[iparams+3], q_reg[i]) 

         circ.cnot(sum(network[:1])+j,i)
      iparams += 4

   return circ,iparams
