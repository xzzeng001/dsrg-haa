from qiskit import QuantumCircuit, QuantumRegister, transpile, Aer, execute
from typing import Union, Optional, List, Tuple, Any, Dict
import numpy as np
from circuit import add_one_qubit_gates, add_two_qubit_gates_1, add_two_qubit_gates_2
from qiskit.providers.aer import QasmSimulator


network=[6,2]
num_qubits=sum(network)
num_params=6*network[0]+4*(network[0]-1)+3*4*(network[0]-1)
params=np.random.uniform(low=-np.pi, high=np.pi, size=(num_params)) #np.loadtxt('params.txt')

q_reg = QuantumRegister(num_qubits, 'q')
# Create a circuit with a register of three qubits
circ = QuantumCircuit(num_qubits)

iparams = 0
circ,iparams = add_one_qubit_gates(circ,q_reg,network,params,iparams)

circ,iparams = add_two_qubit_gates_1(circ,q_reg,network,params,iparams)

circ,iparams = add_two_qubit_gates_2(circ,q_reg,network,params,iparams)

circ,iparams = add_two_qubit_gates_1(circ,q_reg,network,params,iparams)

circ,iparams = add_two_qubit_gates_1(circ,q_reg,network,params,iparams)

circ,iparams = add_one_qubit_gates(circ,q_reg,network,params,iparams)

circ.draw(output='mpl',filename='1.png')

# set the backend, coupling map, basis gates and optimization level of gate_error_probabilities (if given)
transpile_backend = Aer.get_backend('qasm_simulator')
#print('basis gates: ',transpile_backend.configuration().basis_gates)

transpile_coupling_map = None 
transpile_basis_gates = ['h','cz','p']
# optimization level should be ever 0,1,2 or 3
optimization_level = 3

# transpile the quantum circuits
transpiled_circuits = transpile(circ, backend=transpile_backend,
    optimization_level=optimization_level, coupling_map=transpile_coupling_map,
    basis_gates=transpile_basis_gates, seed_transpiler=0)

transpiled_circuits.draw(output='mpl',filename='2.png')
